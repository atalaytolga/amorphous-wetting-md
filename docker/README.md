# HALMD Docker image

The image installs HALMD under `/opt/halmd` and includes Python with `h5py`
for H5MD validation and analysis.

## Build an image

Run builds from the project root. By default, the image is built from the
upstream HALMD repository's `testing` branch.

```bash
docker build \
    --file docker/Dockerfile.halmd \
    --tag halmd \
    docker
```

Pass a repository and ref to build another HALMD revision. For example, build
the NPAT integration branch with:

```bash
docker build \
    --file docker/Dockerfile.halmd \
    --tag halmd:npat-barostat-support \
    --build-arg HALMD_REPOSITORY=https://github.com/atalaytolga/halmd.git \
    --build-arg HALMD_REF=integration/npat-barostat-support \
    docker
```

The Compose file defines four build services. The position-shifting and planar
wall services isolate the two HALMD components for validation; the combined
NPAT service is the production build used by this repository.

| Service | Default image | Default HALMD ref |
| --- | --- | --- |
| `halmd` | `halmd` | `testing` |
| `halmd-position-shifting` | `halmd:position-shifting` | `feature/position-shifting` |
| `halmd-planar-wall-offset-update` | `halmd:planar-wall-offset-update` | `feature/planar-wall-offset-update` |
| `halmd-npat-barostat-support` | `halmd:npat-barostat-support` | `integration/npat-barostat-support` |

Build one service with:

```bash
docker compose -f docker/compose.yaml build halmd-npat-barostat-support
```

The fork builds share `HALMD_FORK_REPOSITORY`. Override their image/ref pairs
with `AWMD_DOCKER_IMAGE_POSITION_SHIFTING` and `HALMD_POSITION_SHIFTING_REF`,
`AWMD_DOCKER_IMAGE_PLANAR_WALL_OFFSET` and `HALMD_PLANAR_WALL_OFFSET_REF`, or
`AWMD_DOCKER_IMAGE_NPAT_BAROSTAT` and `HALMD_NPAT_BAROSTAT_REF`. The upstream
service uses `AWMD_DOCKER_IMAGE`, `HALMD_REPOSITORY`, and `HALMD_REF`.

## Run simulations

Use `scripts/run.sh` for simulations so Git and backend metadata are recorded
in the H5 output:

```bash
./scripts/run.sh \
    --docker-image halmd:npat-barostat-support \
    scripts/simulations/slit_pore/npt.lua \
    --help
```

The backend selectors are:

- `--native` for the host `halmd` executable;
- `--docker` for `halmd`;
- `--docker-position-shifting` for the position-shifting image;
- `--docker-image IMAGE` for any explicit image, including the NPAT image.

Use `--` to stop backend-option parsing when a later argument must be passed
to HALMD.

## Run the container directly

The lower-level Docker wrapper requires an image name. With no command it opens
the image's default shell; otherwise it executes the supplied command:

```bash
AWMD_DOCKER_IMAGE=halmd:npat-barostat-support \
    ./scripts/run/docker/halmd.sh

AWMD_DOCKER_IMAGE=halmd:npat-barostat-support \
    ./scripts/run/docker/halmd.sh halmd --help
```

The wrapper mounts the project at `/work`, forwards `AWMD_*` variables, and
runs with the host user by default. Set `RUN_AS_ROOT=1` to use the container
user, or set `AWMD_DOCKER_GPUS=none` to omit Docker's `--gpus` option.
