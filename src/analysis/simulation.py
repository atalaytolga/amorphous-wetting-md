import h5py
import numpy as np
from typing import Optional, Tuple

class Simulation:
    def __init__(self, path):
        """
        Initialize the Simulation analyzer.
        Args:
            path (str): Path to the HDF5 output file.
        """
        self.path = path
        self._file:  Optional[h5py.File] = None

    def __enter__(self):
        """Context manager entry."""
        try:
            self._file = h5py.File(self.path, 'r')
        except OSError as e:
            raise OSError(f"Could not open file {self.path}. Error: {e}")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Context manager exit (closes file)."""
        _ = (exc_type, exc_value, traceback)

        if self._file:
            self._file.close()
            self._file = None

    def _ensure_open(self):
        """Helper to ensure file is open."""
        if self._file is None:
            raise RuntimeError("File not open. Use 'with Simulation(...) as sim:'")

    @property
    def file(self) -> h5py.File:
        if self._file is None:
            raise RuntimeError("File is not open")
        else:
            return self._file

    def _read_dataset(self, path: str) -> np.ndarray:
        if path not in self.file:
            return np.array([])

        node = self.file[path]

        if isinstance(node, h5py.Dataset):
            return node[()]
        else:
            return np.array([])

    # =========================================================================
    # Geometry & Statistical Helpers
    # =========================================================================

    def _get_box_edges(self) -> np.ndarray:
        self._ensure_open()

        edges = self._read_dataset('particles/all/box/edges')
        return np.diag(edges)

    def _get_box_z(self) -> float:
        edges = self._get_box_edges()
        return edges[2]

    def _compute_block_error(self, data: np.ndarray, n_blocks: int) -> Tuple[np.ndarray, np.ndarray]:
        """Calculates Mean and Standard Error (SEM) using Block Averaging."""
        N = len(data)
        if N == 0:
            return np.array(np.nan), np.array(np.nan)

        mean = np.mean(data, axis=0)

        if n_blocks <= 1 or N < n_blocks:
            std_dev = np.std(data, axis=0, ddof=1)
            return mean, std_dev / np.sqrt(N)

        block_size = N // n_blocks
        cropped_limit = block_size * n_blocks
        cropped_data = data[:cropped_limit]

        shape = (n_blocks, block_size) + cropped_data.shape[1:]
        reshaped = cropped_data.reshape(shape)

        block_means = np.mean(reshaped, axis=1)
        block_std_dev = np.std(block_means, axis=0, ddof=1)
        error = block_std_dev / np.sqrt(n_blocks)

        return mean, error

    def _get_slab_profile(self, observable_name, n_blocks):
        """
        Generic helper to iterate over slab groups in HDF5 and calculate stats.
        """
        self._ensure_open()

        if 'thermo' not in self.file:
            return np.array([]), np.array([])

        thermo = self.file['thermo']

        if not isinstance(thermo, h5py.Group):
            return np.array([]), np.array([])

        # Identify and Sort Slabs Numerically
        slab_keys = [k for k in thermo.keys() if k.startswith('slab')]
        if not slab_keys:
            return np.array([]), np.array([])

        slab_keys.sort(key=lambda x: int(x.split('_')[1]))

        means_list = []
        errors_list = []

        for slab in slab_keys:
            path = f'thermo/{slab}/{observable_name}/value'
            data = self._read_dataset(path)
            if data.size > 0:
                m, e = self._compute_block_error(data, n_blocks)
                means_list.append(m)
                errors_list.append(e)

        return np.array(means_list), np.array(errors_list)

    # =========================================================================
    # Global Observables
    # =========================================================================

    def get_trajectory(self, name):
        """Returns raw time series for a global observable."""
        self._ensure_open()
        base = f'observables/{name}'
        if base not in self.file:
            raise KeyError(f"Observable '{name}' not found.")

        val = self._read_dataset(f'{base}/value')
        time = self._read_dataset(f'{base}/time')
        return val, time

    def get_binned_trajectory(self, name, n_bins=50):
        """Resamples trajectory into n_bins with error bars."""
        try:
            raw_val, raw_time = self.get_trajectory(name)
        except KeyError:
            # Return empty arrays if observable is missing
            return np.array([]), np.array([]), np.array([])

        if len(raw_val) < n_bins:
            return raw_time, raw_val, np.zeros_like(raw_val)

        vals_split = np.array_split(raw_val, n_bins)
        time_split = np.array_split(raw_time, n_bins)

        binned_mean = np.array([np.mean(chunk, axis=0) for chunk in vals_split])
        binned_time = np.array([np.mean(chunk, axis=0) for chunk in time_split])
        binned_std = np.array([np.std(chunk, axis=0) for chunk in vals_split])

        if binned_mean.ndim > 1 and binned_mean.shape[-1] == 1:
            binned_mean = binned_mean.flatten()
            binned_std = binned_std.flatten()

        return binned_time, binned_mean, binned_std

    # =========================================================================
    # Slab Observables (Spatial Profiles)
    # =========================================================================

    def slab_density_profile(self, n_blocks=5):
        """
        Calculate the spatial density profile with Real Z-Coordinates.
        Returns: (z_coords, means, errors)
        """
        # Uses the generic helper to get raw values
        means, errors = self._get_slab_profile('density', n_blocks)

        n_slabs = len(means)
        if n_slabs == 0:
            return np.array([]), np.array([]), np.array([])

        Lz = self._get_box_z()
        dz = Lz / n_slabs
        z_coords = (np.arange(n_slabs) + 0.5) * dz - (Lz / 2.0)

        return z_coords, means, errors

    def slab_stress_profile(self, n_blocks=5):
        """
        Calculates the spatial profile of Diagonal Stress components.
        Returns:
            z_coords (np.array): Center of each slab.
            profiles (dict): Dictionary containing means and errors:
                {
                  'xx': (mean_xx, err_xx),
                  'yy': (mean_yy, err_yy),
                  'zz': (mean_zz, err_zz)
                }
        """
        self._ensure_open()

        Lz = self._get_box_z()

        if 'thermo' not in self.file:
            return None, {}

        thermo = self.file['thermo']

        if not isinstance(thermo, h5py.Group):
             return None, {}

        slab_keys = [k for k in thermo.keys() if k.startswith('slab')]
        if not slab_keys:
            return None, {}

        # Sort numerically
        slab_keys.sort(key=lambda x: int(x.split('_')[1]))

        n_slabs = len(slab_keys)
        dz = Lz / n_slabs
        z_coords = (np.arange(n_slabs) + 0.5) * dz - (Lz / 2.0)

        # 2. Storage
        means = {'xx': [], 'yy': [], 'zz': []}
        errs  = {'xx': [], 'yy': [], 'zz': []}

        # 3. Iterate over Slabs
        for slab in slab_keys:
            path = f'thermo/{slab}/stress_tensor/value'
            data = self._read_dataset(path)
            if data.size > 0:
                data_diag = data[:, 0:3]

                m, e = self._compute_block_error(data_diag, n_blocks)

                # Store (Index 0->xx, 1->yy, 2->zz)
                for i, comp in enumerate(['xx', 'yy', 'zz']):
                    means[comp].append(m[i])
                    errs[comp].append(e[i])

        profiles = {}
        for comp in ['xx', 'yy', 'zz']:
            profiles[comp] = (np.array(means[comp]), np.array(errs[comp]))

        return z_coords, profiles
