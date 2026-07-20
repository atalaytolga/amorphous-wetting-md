"""Create symmetric particle walls from the lower half of an H5MD snapshot."""

import argparse
from pathlib import Path

import h5py
import numpy as np


def get_snapshot(file_handle, group_name="all"):
    """Extract the last position frame and the simulation box."""
    group_path = f"particles/{group_name}"

    if group_path not in file_handle:
        raise KeyError(f"Group '{group_path}' not found.")

    particle_group = file_handle[group_path]
    try:
        position = np.array(particle_group["position"]["value"][-1])
    except KeyError as error:
        raise KeyError(f"Position data missing in '{group_path}'.") from error

    box = np.array(particle_group["box"]["edges"][()])
    return position, box


def main():
    parser = argparse.ArgumentParser(
        description="Symmetrise obstacle walls for HALMD."
    )
    parser.add_argument(
        "--input", required=True, help="Input H5 file containing the initial lattice"
    )
    parser.add_argument("--out", required=True, help="Output H5 file path")
    parser.add_argument("--width", type=float, required=True, help="Pore width")
    parser.add_argument(
        "--group",
        default="obstacles",
        help="Particle group name (default: 'obstacles')",
    )
    args = parser.parse_args()

    if Path(args.out).resolve() == Path(args.input).resolve():
        parser.error("output file must differ from the input file")
    if not np.isfinite(args.width) or args.width <= 0:
        parser.error("width must be finite and positive")

    print(f"Reading input: {args.input}")
    print(f"Target pore width: {args.width}")

    with h5py.File(args.input, "r") as input_file:
        original_position, original_box = get_snapshot(input_file, args.group)

        bottom_position = original_position[original_position[:, 2] < 0.0]
        if len(bottom_position) == 0:
            raise ValueError("No particles found with z < 0; check input coordinates.")

        current_surface = np.max(bottom_position[:, 2])
        target_surface = -args.width / 2.0
        shift_z = target_surface - current_surface

        print(f"   Input surface found at z = {current_surface:.4f}")
        print(
            f"   Shifting bottom wall by dz = {shift_z:.4f} "
            f"to align to z = {target_surface:.4f}"
        )

        bottom_position[:, 2] += shift_z
        top_position = bottom_position.copy()
        top_position[:, 2] *= -1.0
        final_position = np.concatenate((bottom_position, top_position), axis=0)
        final_species = np.ones(len(final_position), dtype=np.uint32)

        print(f"   Final system: {len(final_position)} atoms")
        print(f"   Top surface z    : {np.min(top_position[:, 2]):.4f}")
        print(f"   Bottom surface z : {np.max(bottom_position[:, 2]):.4f}")
        print(f"Writing output to {args.out}...")

        with h5py.File(args.out, "w") as output_file:
            if "h5md" in input_file:
                input_file.copy("h5md", output_file)
            if "parameters" in input_file:
                input_file.copy("parameters", output_file)

            group_path = f"particles/{args.group}"
            particle_group = output_file.create_group(group_path)

            position_group = particle_group.create_group("position")
            position_group.create_dataset(
                "value",
                data=final_position[None, ...].astype(np.float32),
                compression="gzip",
            )
            position_group["step"] = np.array([0], dtype=np.uint64)
            position_group["time"] = np.array([0], dtype=np.float64)

            species_group = particle_group.create_group("species")
            species_group.create_dataset(
                "value",
                data=final_species[None, ...],
                compression="gzip",
            )
            species_group["step"] = np.array([0], dtype=np.uint64)
            species_group["time"] = np.array([0], dtype=np.float64)

            box_group = particle_group.create_group("box")
            box_group.create_dataset("edges", data=original_box.astype(np.float64))

            source_box = input_file[f"{group_path}/box"]
            for attribute in ("boundary", "dimension"):
                if attribute in source_box.attrs:
                    box_group.attrs[attribute] = source_box.attrs[attribute]

    print("Success. Symmetric walls created.")


if __name__ == "__main__":
    main()
