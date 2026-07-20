#!/usr/bin/env bash
set -euo pipefail

# Compare fixed-wall NVT integration with an otherwise identical NPAT run
# that advances the piston on every MD step. Both wall models are active by
# default; set the NPAT_COST_*_WEIGHT variables to benchmark another mixture.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"

source "${PROJECT_ROOT}/scripts/experiments/lib/hdf5.sh"

export AWMD_EXPERIMENT_TYPE="barostat"
export AWMD_EXPERIMENT_NAME="benchmark_npt_cost"

SIMULATION_SCRIPT="scripts/simulations/slit_pore/npt.lua"
DOCKER_IMAGE="${HALMD_IMAGE:-halmd:npat-barostat-support}"
WALL_INPUT="${AWMD_WALL_INPUT:-data/wall/amorphous_symmetrized_W30_rho0.9.h5}"
OUTPUT_DIR="${NPAT_COST_OUTPUT_DIR:-data/benchmarks/npt_cost}"

BENCHMARK_TIME="${NPAT_COST_TIME:-200}"
TIMESTEP="${NPAT_COST_TIMESTEP:-0.005}"
PARTICLE_WEIGHT="${NPAT_COST_PARTICLE_WEIGHT:-0.5}"
FLAT_WEIGHT="${NPAT_COST_FLAT_WEIGHT:-0.5}"

STEPS="$(awk -v time="${BENCHMARK_TIME}" -v dt="${TIMESTEP}" \
    'BEGIN {print int(time / dt + 0.999999999)}')"

particle_walls_enabled() {
    awk -v weight="${PARTICLE_WEIGHT}" 'BEGIN {exit !(weight > 0)}'
}

GEOMETRY_ARGS=(
    --pore_length 20
    --pore_width 30
    --box_width 60
)

SYSTEM_ARGS=(
    --random-seed 20260717
    --fluid_density 0.70
    --temperature 0.8
    --timestep "${TIMESTEP}"
)

SAMPLING_ARGS=(
    --sampling thermodynamics "${STEPS}"
    --sampling barostat "${STEPS}"
    --sampling particles 0
)

WALL_ARGS=(
    --weights particle "${PARTICLE_WEIGHT}"
    --weights flat "${FLAT_WEIGHT}"
)

if particle_walls_enabled; then
    WALL_ARGS=(
        --particle_wall_input "${WALL_INPUT}"
        "${WALL_ARGS[@]}"
    )
fi

BAROSTAT_ARGS=(
    --barostat target_pressure 0.1
    --barostat mass 100000000
    --barostat damping 10
    --barostat max_width_change 0.1
    --barostat interval 1
)

COMMON_ARGS=(
    --overwrite
    "${GEOMETRY_ARGS[@]}"
    "${SYSTEM_ARGS[@]}"
    "${SAMPLING_ARGS[@]}"
    "${WALL_ARGS[@]}"
    "${BAROSTAT_ARGS[@]}"
    --equilibrium_time 0
    --nvt_equilibrium_time 0
)

prepare() {
    if ((STEPS < 1)); then
        printf 'Error: NPAT benchmark requires at least one MD step\n' >&2
        return 2
    fi

    cd -- "${PROJECT_ROOT}"
    require_h5ls
    if particle_walls_enabled; then
        require_readable_h5 "${WALL_INPUT}"
    fi
    mkdir -p -- "${OUTPUT_DIR}"
}

run_case() {
    local label="$1"
    local barostat_time="$2"
    local production_time="$3"
    local result_name="$4"
    local output_basename="${OUTPUT_DIR}/${label}"
    local start_ns
    local end_ns
    local elapsed_seconds

    start_ns="$(date +%s%N)"
    "${PROJECT_ROOT}/scripts/run.sh" \
        --docker-image "${DOCKER_IMAGE}" \
        "${SIMULATION_SCRIPT}" \
        "${COMMON_ARGS[@]}" \
        --output "${output_basename}" \
        --barostat_time "${barostat_time}" \
        --production_time "${production_time}"
    end_ns="$(date +%s%N)"

    elapsed_seconds="$(awk -v start="${start_ns}" -v stop="${end_ns}" \
        'BEGIN {printf "%.9f", (stop - start) / 1e9}')"
    printf -v "${result_name}" '%s' "${elapsed_seconds}"
}

validate_outputs() {
    require_readable_h5 "${OUTPUT_DIR}/fixed_walls.h5"
    require_readable_h5 "${OUTPUT_DIR}/barostat_each_step.h5"
}

md_step_us() {
    awk '/MD integration step:/ {
        for (i = 1; i <= NF; ++i) {
            if ($i == "step:") {
                value = $(i + 1)
                unit = $(i + 2)
            }
        }
    }
    END {
        if (unit == "ns") value /= 1000
        if (unit == "ms") value *= 1000
        if (unit == "s") value *= 1000000
        if (value == "") exit 1
        printf "%.9f", value
    }' "$1"
}

report_results() {
    local fixed_seconds="$1"
    local barostat_seconds="$2"
    local fixed_md_us="$3"
    local barostat_md_us="$4"

    awk \
        -v fixed="${fixed_seconds}" \
        -v barostat="${barostat_seconds}" \
        -v steps="${STEPS}" \
        -v time="${BENCHMARK_TIME}" \
        -v dt="${TIMESTEP}" \
        'BEGIN {
            fixed_us = 1e6 * fixed / steps
            barostat_us = 1e6 * barostat / steps
            slowdown = barostat / fixed
            printf "NPAT cost benchmark: %.6g LJ time, dt=%.6g, %d MD steps\n", time, dt, steps
            printf "fixed walls:           %.3f s total, %.3f us/step\n", fixed, fixed_us
            printf "barostat every step:   %.3f s total, %.3f us/step\n", barostat, barostat_us
            printf "extra cost:            %.3f us/step\n", barostat_us - fixed_us
            printf "slowdown:              %.3fx\n", slowdown
            printf "runtime increase:      %.1f%%\n", 100 * (slowdown - 1)
        }'

    awk \
        -v fixed="${fixed_md_us}" \
        -v barostat="${barostat_md_us}" \
        'BEGIN {
            slowdown = barostat / fixed
            printf "HALMD profiler, fixed:     %.3f us/step\n", fixed
            printf "HALMD profiler, barostat:  %.3f us/step\n", barostat
            printf "HALMD profiler, extra:     %.3f us/step\n", barostat - fixed
            printf "HALMD profiler, slowdown:  %.3fx (%+.1f%%)\n", slowdown, 100 * (slowdown - 1)
        }'
}

main() {
    local fixed_seconds
    local barostat_seconds
    local fixed_md_us
    local barostat_md_us

    prepare

    run_case fixed_walls 0 "${BENCHMARK_TIME}" fixed_seconds
    run_case barostat_each_step "${BENCHMARK_TIME}" 0 barostat_seconds

    validate_outputs

    fixed_md_us="$(md_step_us "${OUTPUT_DIR}/fixed_walls.log")"
    barostat_md_us="$(md_step_us "${OUTPUT_DIR}/barostat_each_step.log")"
    report_results \
        "${fixed_seconds}" \
        "${barostat_seconds}" \
        "${fixed_md_us}" \
        "${barostat_md_us}"
}

main "$@"
