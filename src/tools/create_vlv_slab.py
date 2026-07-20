"""Create a vapour-liquid-vapour slab from two H5MD snapshots."""

import argparse
from pathlib import Path

import h5py
import numpy as np


def wrap_pbc(positions, box_edges):
    """Wrap coordinates into the centred box ``[-L / 2, L / 2]``."""
    box_lengths = np.diag(box_edges)
    return ((positions + box_lengths / 2) % box_lengths) - box_lengths / 2


def remove_com_velocity(velocities):
    """Remove the centre-of-mass velocity to prevent system drift."""
    return velocities - np.mean(velocities, axis=0)


def get_snapshot(file_handle, group_name="all"):
    """Extract the last position and velocity frames and the simulation box."""
    group_path = f"particles/{group_name}"

    if group_path not in file_handle:
        raise KeyError(f"Group '{group_path}' not found in the simulation file.")

    particle_group = file_handle[group_path]

    try:
        position = np.array(particle_group["position"]["value"][-1])
    except KeyError as error:
        raise KeyError(f"Position data missing in '{group_path}'.") from error

    velocity = None
    if "velocity" in particle_group:
        velocity = np.array(particle_group["velocity"]["value"][-1])

    box = np.array(particle_group["box"]["edges"][()])
    return position, velocity, box


def main():
    parser = argparse.ArgumentParser(
        description="Create a vapour-liquid-vapour initial condition for HALMD."
    )
    parser.add_argument(
        "--liq", required=True, help="Input H5 file for the liquid phase"
    )
    parser.add_argument(
        "--vap", required=True, help="Input H5 file for the vapour phase"
    )
    parser.add_argument("--out", required=True, help="Output H5 file path")
    parser.add_argument(
        "--width",
        type=float,
        default=40.0,
        help="Width of the liquid slab in the z direction (default: 40.0)",
    )
    parser.add_argument(
        "--buffer",
        type=float,
        default=0.75,
        help="Safety buffer used to prevent overlap (default: 0.75)",
    )
    parser.add_argument(
        "--group",
        default="all",
        help="Particle group name in the H5 file (default: 'all')",
    )
    args = parser.parse_args()

    output_path = Path(args.out).resolve()
    if output_path in {Path(args.liq).resolve(), Path(args.vap).resolve()}:
        parser.error("output file must differ from both input files")
    if not np.isfinite(args.width) or args.width <= 0:
        parser.error("width must be finite and positive")
    if not np.isfinite(args.buffer) or args.buffer < 0:
        parser.error("buffer must be finite and non-negative")

    print(f"Reading liquid: {args.liq}")
    print(f"Reading vapour: {args.vap}")

    with h5py.File(args.liq, "r") as liquid_file, h5py.File(
        args.vap, "r"
    ) as vapour_file:
        liquid_position, liquid_velocity, liquid_box = get_snapshot(
            liquid_file, args.group
        )
        vapour_position, vapour_velocity, vapour_box = get_snapshot(
            vapour_file, args.group
        )

    has_liquid_velocity = liquid_velocity is not None
    has_vapour_velocity = vapour_velocity is not None
    if has_liquid_velocity != has_vapour_velocity:
        raise ValueError(
            "Velocity data must be present in both input files or absent "
            "from both input files."
        )

    liquid_xy = np.diag(liquid_box)[:2]
    vapour_xy = np.diag(vapour_box)[:2]
    if not np.allclose(vapour_xy, liquid_xy, atol=1e-5):
        raise ValueError(
            "XY dimensions of the liquid and vapour boxes do not match: "
            f"liquid XY = {liquid_xy}, vapour XY = {vapour_xy}."
        )

    z_cut_min = -args.width / 2.0
    z_cut_max = args.width / 2.0

    liquid_position = wrap_pbc(liquid_position, vapour_box)
    vapour_position = wrap_pbc(vapour_position, vapour_box)

    vapour_mask = (vapour_position[:, 2] < z_cut_min - args.buffer) | (
        vapour_position[:, 2] > z_cut_max + args.buffer
    )
    liquid_mask = (liquid_position[:, 2] >= z_cut_min) & (
        liquid_position[:, 2] <= z_cut_max
    )

    selected_liquid_position = liquid_position[liquid_mask]
    selected_vapour_position = vapour_position[vapour_mask]

    print(f"   Selected {len(selected_liquid_position)} liquid atoms")
    print(f"   Selected {len(selected_vapour_position)} vapour atoms")

    combined_position = np.concatenate(
        (selected_liquid_position, selected_vapour_position), axis=0
    )
    combined_velocity = None
    if has_liquid_velocity:
        selected_liquid_velocity = remove_com_velocity(
            liquid_velocity[liquid_mask]
        )
        selected_vapour_velocity = remove_com_velocity(vapour_velocity[vapour_mask])
        combined_velocity = np.concatenate(
            (selected_liquid_velocity, selected_vapour_velocity), axis=0
        )

    print(f"Writing output to {args.out}...")

    with h5py.File(args.out, "w") as output_file:
        with h5py.File(args.liq, "r") as source:
            if "h5md" in source:
                source.copy("h5md", output_file)
            if "parameters" in source:
                source.copy("parameters", output_file)

        group_path = f"particles/{args.group}"
        particle_group = output_file.create_group(group_path)

        position_group = particle_group.create_group("position")
        position_group.create_dataset(
            "value", data=combined_position[None, ...], compression="gzip"
        )
        position_group["step"] = np.array([0], dtype=np.uint64)
        position_group["time"] = np.array([0], dtype=np.float64)

        if combined_velocity is not None:
            velocity_group = particle_group.create_group("velocity")
            velocity_group.create_dataset(
                "value", data=combined_velocity[None, ...], compression="gzip"
            )
            velocity_group["step"] = np.array([0], dtype=np.uint64)
            velocity_group["time"] = np.array([0], dtype=np.float64)

        box_group = particle_group.create_group("box")
        box_group.create_dataset("edges", data=vapour_box)

        with h5py.File(args.liq, "r") as source:
            source_box = source[f"{group_path}/box"]
            for attribute in ("boundary", "dimension"):
                if attribute in source_box.attrs:
                    box_group.attrs[attribute] = source_box.attrs[attribute]

    print("Success")


if __name__ == "__main__":
    main()
