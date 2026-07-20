#!/usr/bin/env bash
set -euo pipefail

# Run fast end-to-end NPAT cases and validate their analysis-oriented H5MD
# output schema.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"

source "${PROJECT_ROOT}/scripts/experiments/lib/hdf5.sh"

export AWMD_EXPERIMENT_TYPE="barostat"
export AWMD_EXPERIMENT_NAME="smoke_npt"

SIMULATION_SCRIPT="scripts/simulations/slit_pore/npt.lua"
VALIDATOR="${PROJECT_ROOT}/src/tests/validate_npt_output.py"
DOCKER_IMAGE="${HALMD_IMAGE:-halmd:npat-barostat-support}"
WALL_INPUT="${AWMD_WALL_INPUT:-data/wall/amorphous_symmetrized_W30_rho0.9.h5}"
OUTPUT_DIR="${NPAT_SMOKE_OUTPUT_DIR:-data/smoke/npt}"
PYTHON="${PYTHON:-python3}"

THERMODYNAMICS_INTERVAL=2
BAROSTAT_SAMPLING_INTERVAL=1
PARTICLE_SAMPLING_INTERVAL=3

prepare() {
    cd -- "${PROJECT_ROOT}"
    require_h5ls
    require_readable_h5 "${WALL_INPUT}"

    if ! command -v "${PYTHON}" >/dev/null 2>&1; then
        printf 'Error: Python interpreter not found: %s\n' "${PYTHON}" >&2
        return 1
    fi
    if [[ ! -f "${VALIDATOR}" ]]; then
        printf 'Error: smoke validator not found: %s\n' "${VALIDATOR}" >&2
        return 1
    fi

    mkdir -p -- "${OUTPUT_DIR}"
}

validate_case() {
    local output_file="$1"
    local label="$2"
    local particle_weight="$3"
    local flat_weight="$4"
    local structure_mode="$5"

    require_readable_h5 "${output_file}"
    "${PYTHON}" "${VALIDATOR}" \
        --input "${output_file}" \
        --label "${label}" \
        --particle-weight "${particle_weight}" \
        --flat-weight "${flat_weight}" \
        --structure-mode "${structure_mode}" \
        --thermodynamics-interval "${THERMODYNAMICS_INTERVAL}" \
        --barostat-sampling-interval "${BAROSTAT_SAMPLING_INTERVAL}" \
        --particle-sampling-interval "${PARTICLE_SAMPLING_INTERVAL}"

    printf 'NPAT smoke passed: %s\n' "${label}"
}

run_case() {
    local label="$1"
    local particle_weight="$2"
    local flat_weight="$3"
    local structure_mode="${4:-none}"
    local output_basename="${OUTPUT_DIR}/${label}"
    local output_file="${output_basename}.h5"
    local fluid_density="0.01"
    local -a analysis_args=()
    local -a simulation_args

    case "${structure_mode}" in
        none)
            ;;
        all|thermo)
            fluid_density="0.05"
            analysis_args=(
                --sampling structure 1
                --sampling structure_observables "${structure_mode}"
                --slab width 7
                --slab axis z
            )
            ;;
        *)
            printf 'Error: unsupported structure mode: %s\n' \
                "${structure_mode}" >&2
            return 2
            ;;
    esac

    simulation_args=(
        # Output and reproducibility
        --output "${output_basename}"
        --overwrite
        --random-seed 20260717

        # Geometry and wall models
        --particle_wall_input "${WALL_INPUT}"
        --pore_length 10
        --pore_width 30
        --box_width 50
        --weights particle "${particle_weight}"
        --weights flat "${flat_weight}"

        # Fluid and integration
        --fluid_density "${fluid_density}"
        --temperature 0.8
        --timestep 0.001

        # Simulation phases
        --equilibrium_time 0.001
        --barostat_time 0.004
        --nvt_equilibrium_time 0.001
        --production_time 0.002

        # Barostat controller
        --barostat target_pressure 0.008
        --barostat mass 1000000
        --barostat damping 10
        --barostat max_width_change 0.001
        --barostat interval 1
        --barostat transition_to_nvt
        --barostat average_time 0.002
        --barostat max_transition_time 0.005
        --barostat switch_tolerance 100

        # Sampling
        --sampling thermodynamics "${THERMODYNAMICS_INTERVAL}"
        --sampling barostat "${BAROSTAT_SAMPLING_INTERVAL}"
        --sampling particles "${PARTICLE_SAMPLING_INTERVAL}"
    )

    "${PROJECT_ROOT}/scripts/run.sh" \
        --docker-image "${DOCKER_IMAGE}" \
        "${SIMULATION_SCRIPT}" \
        "${simulation_args[@]}" \
        "${analysis_args[@]}"

    validate_case \
        "${output_file}" \
        "${label}" \
        "${particle_weight}" \
        "${flat_weight}" \
        "${structure_mode}"
}

main() {
    prepare

    run_case particle_only 1 0
    run_case flat_only 0 1
    run_case particle_and_flat 0.5 0.5
    run_case structure_all 0.5 0.5 all
    run_case structure_thermo 0.5 0.5 thermo
}

main "$@"
