# Amorphous Wetting MD

<img src="docs/halmd-logo.png" alt="HALMD" width="120"/>

Molecular-dynamics simulations for investigating surface behavior of simple fluids on
amorphous solid surfaces, built on HALMD.

## Repository layout

```text
docker/                     HALMD image and Compose builds
docs/                       output-schema documentation
notebooks/                  interactive analysis
scripts/
├── lib/                    reusable Lua modules
├── simulations/            HALMD simulation entrypoints
│   ├── bulk/
│   ├── slit_pore/
│   └── wall_generation/
├── experiments/            executable experiment and smoke-test runners
│   ├── barostat/
│   ├── bulk/
│   └── lib/
└── run.sh                  native/Docker backend dispatcher
src/
├── analysis/               Python analysis utilities
└── tools/                  H5 preparation utilities
```

Lua files under `scripts/simulations` describe individual simulations. Shell
files under `scripts/experiments` provide explicit parameter sets and invoke
those simulations through `scripts/run.sh`.

## Running the slit-pore simulation

The production flat-wall runner uses the NPAT barostat, transitions to a fixed
pore width, and continues with NVT equilibration and production:

```bash
./scripts/experiments/barostat/npt_flat_only.sh
```

The smoke runner checks particle walls, flat walls, combined walls, slab
thermodynamics, and structure-factor output:

```bash
./scripts/experiments/barostat/smoke_npt.sh
```

Both default to the `halmd:npat-barostat-support` image. Set `HALMD_IMAGE` to
override the image. The smoke runner also accepts `AWMD_WALL_INPUT` for the
particle-wall H5 file. A fresh checkout does not contain the ignored
`data/wall/amorphous_symmetrized_W30_rho0.9.h5` default, so provide an existing
wall file through `AWMD_WALL_INPUT` before running the smoke test.

For direct simulation calls, select a backend before the Lua entrypoint:

```bash
./scripts/run.sh \
    --docker-image halmd:npat-barostat-support \
    scripts/simulations/slit_pore/npt.lua \
    --help
```

See [docker/README.md](docker/README.md) for image builds and backend details.

## Output and analysis

Experiment runners require `h5ls` on the host to verify generated H5 files.
Python analysis dependencies can be installed with:

```bash
python3 -m pip install -e .
```

The slit-pore H5 hierarchy is documented in
[docs/h5_schema.md](docs/h5_schema.md). Generated simulation data belongs under
`data/` and is excluded from version control.


## License

Unless otherwise noted, this project is licensed under the GNU General Public
License v3.0 or later (`GPL-3.0-or-later`). See [LICENSE](LICENSE).

Files derived from HALMD retain their original GNU Lesser General Public
License v3.0 or later terms. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for details.