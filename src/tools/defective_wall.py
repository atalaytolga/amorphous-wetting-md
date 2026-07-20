"""Create symmetric hybrid amorphous/crystalline particle walls."""

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


def get_bottom_slab(position, label="input"):
    """Return particles below ``z = 0`` and reject an empty lower slab."""
    bottom_position = position[position[:, 2] < 0.0]
    if len(bottom_position) == 0:
        raise ValueError(
            f"No particles found with z < 0 in the {label} coordinates."
        )
    return bottom_position


def align_slab_surface(position):
    """Shift a slab vertically so that its highest particle is at ``z = 0``."""
    aligned_position = position.copy()
    aligned_position[:, 2] -= np.max(position[:, 2])
    return aligned_position


def main():
    parser = argparse.ArgumentParser(
        description="Create hybrid amorphous/crystalline walls and symmetrize them."
    )
    parser.add_argument(
        "--crystalline",
        required=True,
        help="Input H5 file containing the crystalline base lattice",
    )
    parser.add_argument(
        "--amorphous",
        required=True,
        help="Input H5 file containing the amorphous patch lattice",
    )
    parser.add_argument("--out", required=True, help="Output H5 file path")
    parser.add_argument("--width", type=float, required=True, help="Pore width")
    parser.add_argument(
        "--radius",
        type=float,
        required=True,
        help="Radius of the amorphous patch",
    )
    parser.add_argument(
        "--group",
        default="obstacles",
        help="Particle group name (default: 'obstacles')",
    )
    args = parser.parse_args()

    output_path = Path(args.out).resolve()
    input_paths = {
        Path(args.crystalline).resolve(),
        Path(args.amorphous).resolve(),
    }
    if output_path in input_paths:
        parser.error("output file must differ from both input files")
    if not np.isfinite(args.width) or args.width <= 0:
        parser.error("width must be finite and positive")
    if not np.isfinite(args.radius) or args.radius < 0:
        parser.error("radius must be finite and non-negative")

    print(f"Reading crystalline base: {args.crystalline}")
    print(f"Reading amorphous patch: {args.amorphous}")
    print(f"Target pore width: {args.width}")
    print(f"Patch radius: {args.radius}")

    with h5py.File(args.crystalline, "r") as crystalline_file, h5py.File(
        args.amorphous, "r"
    ) as amorphous_file:
        crystalline_position, crystalline_box = get_snapshot(
            crystalline_file, args.group
        )
        amorphous_position, amorphous_box = get_snapshot(
            amorphous_file, args.group
        )

        crystalline_xy = np.diag(crystalline_box)[:2]
        amorphous_xy = np.diag(amorphous_box)[:2]
        if not np.allclose(crystalline_xy, amorphous_xy, atol=1e-5):
            raise ValueError(
                "XY dimensions of the crystalline and amorphous boxes do not "
                "match: "
                f"crystalline XY = {crystalline_xy}, "
                f"amorphous XY = {amorphous_xy}."
            )

        crystalline_bottom = align_slab_surface(
            get_bottom_slab(crystalline_position, "crystalline")
        )
        amorphous_bottom = align_slab_surface(
            get_bottom_slab(amorphous_position, "amorphous")
        )
        print("   Aligned input slabs to z_surface = 0.0 for merging")

        crystalline_radius_squared = (
            crystalline_bottom[:, 0] ** 2 + crystalline_bottom[:, 1] ** 2
        )
        amorphous_radius_squared = (
            amorphous_bottom[:, 0] ** 2 + amorphous_bottom[:, 1] ** 2
        )
        patch_radius_squared = args.radius**2

        keep_crystalline = crystalline_radius_squared >= patch_radius_squared
        keep_amorphous = amorphous_radius_squared < patch_radius_squared
        hybrid_bottom = np.concatenate(
            (
                crystalline_bottom[keep_crystalline],
                amorphous_bottom[keep_amorphous],
            ),
            axis=0,
        )
        if len(hybrid_bottom) == 0:
            raise ValueError("The selected amorphous/crystalline wall is empty.")

        print(
            f"   Merged: {np.sum(keep_amorphous)} amorphous atoms + "
            f"{np.sum(keep_crystalline)} crystalline atoms"
        )

        current_surface = np.max(hybrid_bottom[:, 2])
        target_surface = -args.width / 2.0
        shift_z = target_surface - current_surface

        print(f"   Hybrid surface currently at z = {current_surface:.4f}")
        print(
            f"   Shifting bottom wall by dz = {shift_z:.4f} "
            f"to align to z = {target_surface:.4f}"
        )

        hybrid_bottom[:, 2] += shift_z
        top_position = hybrid_bottom.copy()
        top_position[:, 2] *= -1.0

        final_position = np.concatenate((hybrid_bottom, top_position), axis=0)
        final_species = np.ones(len(final_position), dtype=np.uint32)

        print(f"   Final system: {len(final_position)} atoms")
        print(f"   Top surface z    : {np.min(top_position[:, 2]):.4f}")
        print(f"   Bottom surface z : {np.max(hybrid_bottom[:, 2]):.4f}")
        print(f"Writing output to {args.out}...")

        with h5py.File(args.out, "w") as output_file:
            if "h5md" in crystalline_file:
                crystalline_file.copy("h5md", output_file)
            if "parameters" in crystalline_file:
                crystalline_file.copy("parameters", output_file)

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
            box_group.create_dataset(
                "edges", data=crystalline_box.astype(np.float64)
            )

            source_box = crystalline_file[f"{group_path}/box"]
            for attribute in ("boundary", "dimension"):
                if attribute in source_box.attrs:
                    box_group.attrs[attribute] = source_box.attrs[attribute]

    print("Success. Hybrid symmetric walls created.")


if __name__ == "__main__":
    main()
