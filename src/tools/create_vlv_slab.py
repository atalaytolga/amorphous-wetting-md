import h5py
import numpy as np
import argparse
import sys

def wrap_pbc(positions, box_edges):
    """
    Wraps coordinates into the centered box [-L/2, L/2].
    """
    # Extract diagonal L from box edges
    L = np.diag(box_edges)

    return ((positions + L/2) % L) - L/2

def remove_com_velocity(velocities):
     """
     Remove center-off-mass velocities to prevent system drift.
     """
     return velocities - np.mean(velocities, axis = 0)


def get_snapshot(file_handle, group_name="all"):
    """
    Extracts the last frame of position, velocity, and box from H5MD file.
    """
    group_path = f"particles/{group_name}"

    if group_path not in file_handle:
        raise KeyError(f"Group '{group_path} not found in the simulation file.")

    p_group = file_handle[f"particles/{group_name}"]

    try:
        pos = np.array(p_group["position"]["value"][-1])
    except KeyError:
        raise KeyError(f"Position data missing in {group_path}")

    vel = None
    if "velocity" in p_group:
        vel = np.array(p_group["velocity"]["value"][-1])

    box = np.array(p_group["box"]["edges"][()])

    return pos, vel, box

def main():
    parser = argparse.ArgumentParser(description="Create a Vapor-Liquid-Vapor initial condition for HALMD.")

    # Required Arguments
    parser.add_argument("--liq", required=True, help="Input H5 file for LIQUID phase")
    parser.add_argument("--vap", required=True, help="Input H5 file for VAPOR phase")
    parser.add_argument("--out", required=True, help="Output H5 file path")

    # Optional Arguments
    parser.add_argument("--width", type=float, default=40.0, help="Width of the liquid slab in Z direction (default: 40.0)")
    parser.add_argument("--buffer", type=float, default=0.75, help="Safety buffer to prevent overlap (default: 0.75)")
    parser.add_argument("--group", type=str, default="all", help="Particle group name in H5 file (default: 'all')")

    args = parser.parse_args()

    print(f"Reading Liquid: {args.liq}")
    print(f"Reading Vapor:  {args.vap}")

    with h5py.File(args.liq, "r") as f_liq, h5py.File(args.vap, 'r') as f_vap:
            pos_liq, vel_liq, box_liq = get_snapshot(f_liq, args.group)
            pos_vap, vel_vap, box_vap = get_snapshot(f_vap, args.group)


    box_diff = np.diag(box_vap)[:2] - np.diag(box_liq)[:2]

    if not np.allclose(box_diff, 0, atol=1e-5):
            liq_dims = np.diag(box_liq[:2])
            vap_dims = np.diag(box_vap[:2])
            raise ValueError(
               f"Error: XY dimensions of Liquid and Vapor boxes do not match!"
               f"Liquid XY: {liq_dims}"
               f"Vapor XY:  {vap_dims}"
            )


    z_cut_min = - (args.width / 2.0)
    z_cut_max = + (args.width / 2.0)

    pos_liq = wrap_pbc(pos_liq, box_vap)
    pos_vap = wrap_pbc(pos_vap, box_vap)

    z_vap_min = z_cut_min - args.buffer
    z_vap_max = z_cut_max + args.buffer

    mask_vap = (pos_vap[:, 2] < z_vap_min) | (pos_vap[:, 2] > z_vap_max)
    mask_liq = (pos_liq[:, 2] >= z_cut_min) & (pos_liq[:, 2] <= z_cut_max)

    final_pos_vap = pos_vap[mask_vap]
    final_pos_liq = pos_liq[mask_liq]

    print(f"   Selected {len(final_pos_liq)} Liquid atoms")
    print(f"   Selected {len(final_pos_vap)} Vapor atoms")

    # Concatenate Positions
    combined_pos = np.concatenate((final_pos_liq, final_pos_vap), axis=0)
    combined_vel = None

    if vel_vap is not None and vel_liq is not None:
        final_vel_liq = vel_liq[mask_liq]
        final_vel_vap = vel_vap[mask_vap]

        # Remove COM velocities:
        final_vel_liq = remove_com_velocity(final_vel_liq)
        final_vel_vap = remove_com_velocity(final_vel_vap)

        # Concatenate Velocities
        combined_vel = np.concatenate((final_vel_liq, final_vel_vap), axis=0)

    # Create Species Tags (Change np.zeros to np.ones if two species are required)
    species_liq = np.zeros(len(final_pos_liq), dtype=np.int32)
    species_vap = np.zeros(len(final_pos_vap), dtype=np.int32)
    combined_species = np.concatenate((species_liq, species_vap), axis=0)

    print(f"Writing output to {args.out}...")

    with h5py.File(args.out, 'w') as f_out:
        # Copy H5MD structural metadata from one of the source files
        with h5py.File(args.liq, "r") as source:
            if "h5md" in source:
                source.copy("h5md", f_out)
            if "parameters" in source:
                source.copy("parameters", f_out)

        #if "parameters" in source.keys():
        #    source.copy("parameters", f_out)

        # Create particles group
        group_path = f"particles/{args.group}"
        p_all = f_out.create_group(group_path)

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

        # Copy Box attributes (dimension, boundary conditions)
        with h5py.File(args.liq, "r") as src:
            src_box = src[f"{group_path}/box"]
            for attr in ["boundary", "dimension"]:
                if attr in src_box.attrs:
                    p_box.attrs[attr] = src_box.attrs[attr]
    print("Success")

if __name__ == "__main__":
    main()
