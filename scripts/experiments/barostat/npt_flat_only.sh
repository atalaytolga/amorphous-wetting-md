#!/usr/bin/env bash
set -euo pipefail

# Run a production NPAT-to-NVT simulation confined only by flat walls.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"

source "${PROJECT_ROOT}/scripts/experiments/lib/hdf5.sh"

export AWMD_EXPERIMENT_TYPE="barostat"
export AWMD_EXPERIMENT_NAME="npt_flat_only"

SIMULATION_SCRIPT="scripts/simulations/slit_pore/npt.lua"
DOCKER_IMAGE="${HALMD_IMAGE:-halmd:npat-barostat-support}"
OUTPUT_DIR="${NPAT_FLAT_OUTPUT_DIR:-data/production}"
OUTPUT_NAME="${NPAT_FLAT_OUTPUT_NAME:-npat_flat_run1}"
OUTPUT_BASENAME="${OUTPUT_DIR}/${OUTPUT_NAME}"
OUTPUT_FILE="${OUTPUT_BASENAME}.h5"

prepare() {
    cd -- "${PROJECT_ROOT}"
    require_h5ls
    mkdir -p -- "${OUTPUT_DIR}"
}

run_simulation() {
    local -a simulation_args=(
        # Output and reproducibility
        --output "${OUTPUT_BASENAME}"
        --overwrite
        --random-seed 20260718

        # Geometry and wall models
        --pore_length 20
        --pore_width 60
        --box_width 80
        --weights particle 0
        --weights flat 1

        # Fluid and integration
        --fluid_density 0.70
        --temperature 0.8
        --timestep 0.005

        # Simulation phases
        --equilibrium_time 250
        --barostat_time 500
        --nvt_equilibrium_time 250
        --production_time 500

        # Barostat controller
        --barostat target_pressure 0.1
        --barostat mass 3200
        --barostat damping 800
        --barostat max_width_change 0.1
        --barostat interval 1
        --barostat transition_to_nvt
        --barostat average_time 10
        --barostat max_transition_time 500
        --barostat switch_tolerance 0.02

        # Sampling and analysis
        --sampling thermodynamics 1000
        --sampling barostat 1000
        --sampling particles 1000
        --sampling structure 1000
        --sampling structure_observables all
        --slab width 0.20
        --slab axis z
    )

    "${PROJECT_ROOT}/scripts/run.sh" \
        --docker-image "${DOCKER_IMAGE}" \
        "${SIMULATION_SCRIPT}" \
        "${simulation_args[@]}"
}

validate_output() {
    require_readable_h5 "${OUTPUT_FILE}"
}

main() {
    prepare
    run_simulation
    validate_output
    printf 'NPAT flat-wall output: %s\n' "${OUTPUT_FILE}"
}

main "$@"
