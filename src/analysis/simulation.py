from typing import Optional, Tuple

import h5py
import numpy as np


class Simulation:
    """Read and analyse one HALMD HDF5 output file."""

    def __init__(self, path):
        """Initialize the analyzer for ``path``."""
        self.path = path
        self._file: Optional[h5py.File] = None

    def __enter__(self):
        """Open the HDF5 file for reading."""
        try:
            self._file = h5py.File(self.path, "r")
        except OSError as error:
            raise OSError(f"Could not open {self.path}: {error}") from error
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Close the HDF5 file."""
        _ = (exc_type, exc_value, traceback)
        if self._file is not None:
            self._file.close()
            self._file = None

    def _ensure_open(self):
        if self._file is None:
            raise RuntimeError("File not open. Use 'with Simulation(...) as sim:'")

    @property
    def file(self) -> h5py.File:
        """Return the open HDF5 file handle."""
        self._ensure_open()
        assert self._file is not None
        return self._file

    def _read_dataset(self, path: str) -> np.ndarray:
        path = path.lstrip("/")
        if path not in self.file:
            return np.array([])

        node = self.file[path]
        if isinstance(node, h5py.Dataset):
            return np.asarray(node[()])

        if isinstance(node, h5py.Group) and "value" in node:
            value = node["value"]
            if isinstance(value, h5py.Dataset):
                return np.asarray(value[()])

        return np.array([])

    def read_observable(self, name: str, skip_frames: int = 0) -> np.ndarray:
        """Read an element below ``/observables``.

        ``name`` may include nested groups, for example
        ``"global/thermodynamics/pressure"``.
        """
        if skip_frames < 0:
            raise ValueError("skip_frames must be non-negative")

        data = self._read_dataset(f"observables/{name}")
        if data.size == 0:
            raise KeyError(f"Observable {name!r} not found")

        if data.ndim == 0:
            if skip_frames:
                raise ValueError("Cannot skip frames from a scalar observable")
            return data

        return data[skip_frames:]

    def read_observable_scalar(self, name: str) -> float:
        """Read an observable stored as one scalar value."""
        data = np.squeeze(self.read_observable(name))
        if data.ndim != 0:
            raise ValueError(
                f"Observable {name!r} is not scalar; found shape {data.shape}"
            )
        return float(data)

    def observable_block_means(
        self,
        name: str,
        n_blocks: int,
        discard_fraction: float = 0.0,
    ) -> np.ndarray:
        """Return non-overlapping block means for a time-series observable."""
        if n_blocks < 1:
            raise ValueError("n_blocks must be positive")
        if not 0 <= discard_fraction < 1:
            raise ValueError("discard_fraction must be in [0, 1)")

        data = np.squeeze(self.read_observable(name))
        if data.ndim != 1:
            raise ValueError(
                f"Observable {name!r} must be one-dimensional; "
                f"found shape {data.shape}"
            )

        retained = data[int(discard_fraction * data.size) :]
        block_size = retained.size // n_blocks
        if block_size < 1:
            raise ValueError(
                f"Only {retained.size} retained samples for {n_blocks} blocks"
            )

        retained = retained[: n_blocks * block_size]
        return retained.reshape(n_blocks, block_size).mean(axis=1)

    # Geometry and statistical helpers

    @staticmethod
    def _validate_axis(axis: str) -> str:
        if axis not in {"x", "y", "z"}:
            raise ValueError("axis must be 'x', 'y', or 'z'")
        return axis

    @staticmethod
    def _validate_sample_selection(
        skip_frames: int,
        discard_fraction: float,
    ) -> None:
        if skip_frames < 0:
            raise ValueError("skip_frames must be non-negative")
        if not 0 <= discard_fraction < 1:
            raise ValueError("discard_fraction must be in [0, 1)")
        if skip_frames and discard_fraction:
            raise ValueError("Use either skip_frames or discard_fraction, not both")

    def _get_box_edges(self, group=None) -> np.ndarray:
        """Return the box side lengths from canonical or particle metadata."""
        paths = ["geometry/box/edges"]
        if group is not None:
            paths.append(f"particles/{group}/box/edges")

        for path in paths:
            edges = self._read_dataset(path)
            if edges.shape == (3, 3):
                return np.diag(edges)
        return np.array([])

    def _get_box_length(self, group, axis: str, fallback=None) -> float:
        axis = self._validate_axis(axis)
        edges = self._get_box_edges(group)
        axis_index = {"x": 0, "y": 1, "z": 2}[axis]
        if len(edges) > axis_index:
            return float(edges[axis_index])
        if fallback is not None:
            return float(fallback)
        raise ValueError(
            f"Could not determine the box {axis}-length for group={group!r}; "
            "box edges are missing from the HDF5 file."
        )

    def _get_box_z(self, group, fallback=None) -> float:
        """Compatibility wrapper for callers using the former z-only helper."""
        return self._get_box_length(group, "z", fallback)

    def _get_slab_group(self, axis: str) -> Optional[h5py.Group]:
        axis = self._validate_axis(axis)
        path = f"observables/slabs/{axis}"
        if path not in self.file:
            return None
        group = self.file[path]
        return group if isinstance(group, h5py.Group) else None

    def _get_slab_keys(self, axis: str) -> list:
        group = self._get_slab_group(axis)
        if group is None:
            return []

        keys = [
            key
            for key, node in group.items()
            if isinstance(node, h5py.Group)
            and key.startswith("slab_")
            and key[5:].isdigit()
        ]
        return sorted(keys, key=lambda key: int(key[5:]))

    def _get_slab_centers(
        self,
        axis: str,
        n_slabs: int,
        group=None,
        fallback=None,
    ) -> np.ndarray:
        """Read slab centers from schema-v2 geometry metadata.

        The box-based calculation is retained only as a compatibility fallback
        when ``/geometry/slabs/<axis>`` metadata is absent.
        """
        axis = self._validate_axis(axis)
        geometry_path = f"geometry/slabs/{axis}"
        if geometry_path in self.file:
            geometry = self.file[geometry_path]
            if isinstance(geometry, h5py.Group) and "center" in geometry.attrs:
                centers = np.atleast_1d(
                    np.asarray(geometry.attrs["center"], dtype=float)
                )
                declared_count = geometry.attrs.get("count")
                if declared_count is not None and int(declared_count) != len(centers):
                    raise ValueError(
                        f"{geometry_path} declares {int(declared_count)} slabs "
                        f"but stores {len(centers)} centers"
                    )
                if len(centers) != n_slabs:
                    raise ValueError(
                        f"{geometry_path} stores {len(centers)} centers but "
                        f"observables contain {n_slabs} slabs"
                    )
                return centers

        length = self._get_box_length(group, axis, fallback)
        slab_width = length / n_slabs
        return (np.arange(n_slabs) + 0.5) * slab_width - 0.5 * length

    def _get_slab_series(
        self,
        observable_name: str,
        axis: str = "z",
        skip_frames: int = 0,
        discard_fraction: float = 0.0,
    ) -> np.ndarray:
        """Read one thermodynamic time series from every slab."""
        axis = self._validate_axis(axis)
        self._validate_sample_selection(skip_frames, discard_fraction)
        slab_keys = self._get_slab_keys(axis)
        if not slab_keys:
            return np.array([])

        series = []
        for slab in slab_keys:
            path = (
                f"observables/slabs/{axis}/{slab}/thermodynamics/"
                f"{observable_name}"
            )
            data = self._read_dataset(path)
            if data.size == 0:
                return np.array([])
            if data.ndim == 0:
                raise ValueError(f"Slab observable {path!r} is not time-dependent")

            start = (
                int(discard_fraction * len(data))
                if discard_fraction
                else skip_frames
            )
            series.append(data[start:])

        sample_counts = {len(data) for data in series}
        if len(sample_counts) != 1:
            raise ValueError(
                f"Slab observable {observable_name!r} has inconsistent sample counts"
            )
        if not sample_counts or next(iter(sample_counts)) == 0:
            return np.array([])
        return np.stack(series)

    def _compute_block_error(
        self,
        data: np.ndarray,
        n_blocks: int,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate a mean and standard error using block averaging."""
        if n_blocks < 1:
            raise ValueError("n_blocks must be positive")

        sample_count = len(data)
        if sample_count == 0:
            return np.array(np.nan), np.array(np.nan)

        mean = np.mean(data, axis=0)
        if sample_count == 1:
            return mean, np.zeros_like(mean, dtype=float)

        if n_blocks == 1 or sample_count < n_blocks:
            standard_deviation = np.std(data, axis=0, ddof=1)
            return mean, standard_deviation / np.sqrt(sample_count)

        block_size = sample_count // n_blocks
        cropped_data = data[: block_size * n_blocks]
        shape = (n_blocks, block_size) + cropped_data.shape[1:]
        block_means = np.mean(cropped_data.reshape(shape), axis=1)
        block_standard_deviation = np.std(block_means, axis=0, ddof=1)
        return mean, block_standard_deviation / np.sqrt(n_blocks)

    def _get_slab_profile(
        self,
        observable_name,
        n_blocks,
        skip_frames=0,
        discard_fraction=0.0,
        axis="z",
    ):
        """Return per-slab means and block-averaged standard errors."""
        series = self._get_slab_series(
            observable_name,
            axis=axis,
            skip_frames=skip_frames,
            discard_fraction=discard_fraction,
        )
        if series.size == 0:
            return np.array([]), np.array([])

        statistics = [
            self._compute_block_error(slab_series, n_blocks)
            for slab_series in series
        ]
        means, errors = zip(*statistics)
        return np.asarray(means), np.asarray(errors)

    # Global observables

    def get_trajectory(self, name, skip_frames=0):
        """Return a global observable and its time axis."""
        base = f"observables/{name}"
        if base not in self.file:
            raise KeyError(f"Observable {name!r} not found")

        values = self.read_observable(name, skip_frames=skip_frames)
        time = self._read_dataset(f"{base}/time")
        if skip_frames:
            time = time[skip_frames:]
        return values, time

    def get_binned_trajectory(self, name, n_bins=50, skip_frames=0):
        """Resample a trajectory into equally sized consecutive bins."""
        if n_bins < 1:
            raise ValueError("n_bins must be positive")

        try:
            raw_values, raw_time = self.get_trajectory(
                name,
                skip_frames=skip_frames,
            )
        except KeyError:
            return np.array([]), np.array([]), np.array([])

        sample_count = len(raw_values)
        if sample_count < n_bins:
            return raw_time, raw_values, np.zeros_like(raw_values)

        block_size = sample_count // n_bins
        limit = block_size * n_bins
        value_shape = (n_bins, block_size) + raw_values.shape[1:]
        time_shape = (n_bins, block_size)
        reshaped_values = raw_values[:limit].reshape(value_shape)
        reshaped_time = raw_time[:limit].reshape(time_shape)

        binned_mean = np.mean(reshaped_values, axis=1)
        binned_time = np.mean(reshaped_time, axis=1)
        if block_size > 1:
            binned_standard_deviation = np.std(
                reshaped_values,
                axis=1,
                ddof=1,
            )
        else:
            binned_standard_deviation = np.zeros_like(binned_mean)

        if binned_mean.ndim > 1 and binned_mean.shape[-1] == 1:
            binned_mean = binned_mean.flatten()
            binned_standard_deviation = binned_standard_deviation.flatten()

        return binned_time, binned_mean, binned_standard_deviation

    # Slab observables

    def slab_density_profile(
        self,
        group=None,
        n_blocks=5,
        skip_frames=0,
        fallback=None,
        discard_fraction=0.0,
        axis="z",
    ):
        """Return slab centers, mean densities, and their standard errors.

        ``group`` and ``fallback`` are retained for caller compatibility.
        Schema-v2 files use ``/geometry/slabs/<axis>@center``.
        """
        means, errors = self._get_slab_profile(
            "density",
            n_blocks,
            skip_frames=skip_frames,
            discard_fraction=discard_fraction,
            axis=axis,
        )
        if means.size == 0:
            return np.array([]), np.array([]), np.array([])

        centers = self._get_slab_centers(
            axis,
            len(means),
            group=group,
            fallback=fallback,
        )
        return centers, means, errors

    def slab_stress_profile(
        self,
        group=None,
        n_blocks=5,
        skip_frames=0,
        fallback=None,
        discard_fraction=0.0,
        axis="z",
    ):
        """Return slab centers and profiles of the diagonal stress components."""
        means, errors = self._get_slab_profile(
            "stress_tensor",
            n_blocks,
            skip_frames=skip_frames,
            discard_fraction=discard_fraction,
            axis=axis,
        )
        if means.size == 0:
            return None, {}
        if means.ndim != 2 or means.shape[1] < 3:
            raise ValueError(
                "Slab stress tensors must contain xx, yy, and zz components"
            )

        centers = self._get_slab_centers(
            axis,
            len(means),
            group=group,
            fallback=fallback,
        )
        profiles = {
            component: (means[:, index], errors[:, index])
            for index, component in enumerate(("xx", "yy", "zz"))
        }
        return centers, profiles

    def liquid_vapour_surface_tension(
        self,
        group,
        area,
        n_blocks=5,
        skip_frames=0,
        discard_fraction=0.0,
        axis="z",
    ):
        """Calculate surface tension from slab-resolved stress tensors.

        The leading ``group`` argument is retained for API compatibility and is
        not needed for schema-v2 output.
        """
        _ = group
        axis = self._validate_axis(axis)
        if not np.isfinite(area) or area <= 0:
            raise ValueError("area must be finite and positive")

        stress = self._get_slab_series(
            "stress_tensor",
            axis=axis,
            skip_frames=skip_frames,
            discard_fraction=discard_fraction,
        )
        if stress.size == 0:
            raise KeyError(f"No slab stress tensors found for axis {axis!r}")
        if stress.ndim != 3 or stress.shape[2] < 3:
            raise ValueError(
                "Slab stress tensors must contain xx, yy, and zz components"
            )

        normal_index = {"x": 0, "y": 1, "z": 2}[axis]
        tangential_indices = [index for index in range(3) if index != normal_index]
        normal_stress = stress[:, :, normal_index]
        tangential_stress = np.mean(stress[:, :, tangential_indices], axis=2)
        gamma_time_series = np.sum(normal_stress - tangential_stress, axis=0)
        gamma_time_series /= 2 * area

        gamma_mean, gamma_error = self._compute_block_error(
            gamma_time_series,
            n_blocks,
        )
        return gamma_mean, gamma_error, gamma_time_series
