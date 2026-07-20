#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
SIMULATION_SCRIPT="scripts/simulations/bulk/homogeneous.lua"

export AWMD_EXPERIMENT_TYPE="bulk"
export AWMD_EXPERIMENT_NAME="homogeneous"

main() {
    exec "${PROJECT_ROOT}/scripts/run.sh" "${SIMULATION_SCRIPT}" "$@"
}

main "$@"
