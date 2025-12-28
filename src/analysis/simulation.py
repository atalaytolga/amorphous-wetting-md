import h5py
import numpy as np

class Simulation:
    def __init__(self, path):
        self.path = path
        self._file = None

    def __enter__(self):
        self._file = h5py.File(self.path, 'r')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._file:
            self._file.close()

    def potential_energy(self):
        val = np.array(self._file['observables/potential_energy/value'])
        time = np.array(self._file['observables/potential_energy/time'])
        return np.stack((val, time))
    
    def pressure(self):
        val = np.array(self._file['observables/pressure/value'])
        time = np.array(self._file['observables/pressure/time'])
        return np.stack((val, time))
    
    def temperature(self):
        val = np.array(self._file['observables/temperature/value'])
        time = np.array(self._file['observables/temperature/time'])
        return np.stack((val, time))

    def stress_tensor(self):
        return np.array(self._file['observables/stress_tensor/value'])
    
    def slab_stress_tensor_sampled(self):
        """Calculate time averaged stress tensor over slabs."""
        # Get list of slabs (e.g. slab_1, slab_2...)
        region_keys = [k for k in self._file['thermo'].keys() if k.startswith('slab')]
        num_regions = len(region_keys)
    
        # Initialize: (Number of slabs, 6 tensor components)
        stress_tensor_sampling = np.zeros((num_regions, 6))
        
        for i in range(1, num_regions + 1):
            dataset_path = f'thermo/slab_{i}/stress_tensor/value'
            if dataset_path in self._file:
                cur_data = self._file[dataset_path]
                # CRITICAL FIX: axis=0 averages over time (rows), keeping the 6 components
                stress_tensor_sampling[i - 1] = np.mean(cur_data, axis=0)
        
        return stress_tensor_sampling

    def slab_density_sampled(self):
        """Calculate time averaged density over slabs."""
        region_keys = [k for k in self._file['thermo'].keys() if k.startswith('slab')]
        num_regions = len(region_keys)
    
        density_sampling = np.zeros(num_regions)
        
        for i in range(1, num_regions + 1):
            dataset_path = f'thermo/slab_{i}/density/value'
            if dataset_path in self._file:
                cur_data = self._file[dataset_path]
                density_sampling[i - 1] = np.mean(cur_data, axis=0)
        
        return density_sampling