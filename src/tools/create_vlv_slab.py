import h5py
import numpy as np

# --- CONFIGURATION ---
FILE_LIQ = "#"
FILE_VAP = "#"
FILE_OUT = "#"

GAP_Z_WIDTH = 40.0
BUFFER = 0.75  # Safety distance to prevent overlap at the interface


def wrap_pbc(positions, box_edges):
    """
    Wraps coordinates into the centered box [-L/2, L/2].
    """
    # Extract diagonal L from box edges
    L = np.diag(box_edges)

    return ((positions + L/2) % L) - L/2

def remove_com_velocity(velocities):
     return velocities - np.mean(velocities, axis = 0)


def get_snapshot(file_handle):
    """
    Extracts the last frame of position, velocity and the box from H5MD file.
    """

    p_group = file_handle["particles"]["all"] # Adjust group name if not 'all'

    # Get the last frame [-1]

    pos = p_group["position"]["value"][-1]
    vel = p_group["velocity"]["value"][-1]
    box = p_group["box"]["edges"][()]

    return pos, vel, box

print("Loading Data...")
with h5py.File(FILE_LIQ, "r") as f_liq, h5py.File(FILE_VAP, 'r') as f_vap:
    pos_liq, vel_liq, box_liq = get_snapshot(f_liq)
    pos_vap, vel_vap, box_vap = get_snapshot(f_vap)

assert np.all(np.diag(box_vap)[:2] == np.diag(box_liq)[:2]), "XY dimensions mismatch"

z_cut_min = - (GAP_Z_WIDTH / 2.0)
z_cut_max = + (GAP_Z_WIDTH / 2.0)

# Use Vapor box as Master
L_master = np.diag(box_vap)

pos_liq = wrap_pbc(pos_liq, box_vap)
pos_vap = wrap_pbc(pos_vap, box_vap)

z_vap_min = z_cut_min - BUFFER
z_vap_max = z_cut_max + BUFFER

mask_vap = (pos_vap[:, 2] < z_vap_min) | (pos_vap[:, 2] > z_vap_max)
final_pos_vap = pos_vap[mask_vap]
final_vel_vap = vel_vap[mask_vap]


mask_liq = (pos_liq[:, 2] >= z_cut_min) & (pos_liq[:, 2] <= z_cut_max)
final_pos_liq = pos_liq[mask_liq]
final_vel_liq = vel_liq[mask_liq]


print(f"   Selected {len(final_pos_liq)} Liquid atoms")
print(f"   Selected {len(final_pos_vap)} Vapor atoms")


# Concatenate Positions
combined_pos = np.concatenate((final_pos_liq, final_pos_vap), axis=0)

# Remove COM velocities:
final_vel_liq = remove_com_velocity(final_vel_liq)
final_vel_vap = remove_com_velocity(final_vel_vap)

# Concatenate Velocities
combined_vel = np.concatenate((final_vel_liq, final_vel_vap), axis=0)

# Create Species Tags (Change np.zeros to np.ones if two species are required)
species_liq = np.zeros(len(final_pos_liq), dtype=np.int32)
species_vap = np.zeros(len(final_pos_vap), dtype=np.int32)
combined_species = np.concatenate((species_liq, species_vap), axis=0)

print(f"Writing {FILE_OUT}...")
with h5py.File(FILE_OUT, 'w') as f_out:

    # Copy H5MD Header from Liquid Simulation
    with h5py.File(FILE_LIQ, "r") as source:
        source.copy("h5md", f_out)

    #if "parameters" in source.keys():
    #    source.copy("parameters", f_out)

    # Create particles group
    p_all = f_out.create_group("particles/all")

    # Copy Positions (Time-dependent shape: [1, N, 3])
    p_pos = p_all.create_group("position")
    p_pos.create_dataset("value", data=combined_pos[None, ...], compression="gzip")
    p_pos["step"] = np.array([0], dtype=np.int32)
    p_pos["time"] = np.array([0], dtype=np.int32)

    # Copy Velocities (Time-dependent shape: [1, N, 3])
    p_vel = p_all.create_group("velocity")
    p_vel.create_dataset("value", data=combined_vel[None, ...], compression="gzip")
    p_vel["step"] = np.array([0], dtype=np.int32)
    p_vel["time"] = np.array([0], dtype=np.int32)

    # Create species 
    # TO-DO

    # Create Box
    p_box = p_all.create_group("box")
    p_box.create_dataset("edges", data=box_vap)

    with h5py.File(FILE_LIQ, "r") as src:
        src_box = src["particles/all/box"]
        if "boundary" in src_box.attrs:
            p_box.attrs["boundary"] = src_box.attrs["boundary"]
        if "dimension" in src_box.attrs:
            p_box.attrs["dimension"] = src_box.attrs["dimension"]
print("Done")