#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

source "${PROJECT_ROOT}/scripts/run/lib/backend.sh"
source "${PROJECT_ROOT}/scripts/run/lib/metadata_env.sh"

main() {
    parse_backend_args "$@"

    if (( ${#BACKEND_ARGS[@]} == 0 )); then
        printf 'Error: no simulation script specified\n' >&2
        return 2
    fi

    export AWMD_SIMULATION_SCRIPT="${BACKEND_ARGS[0]}"

    set_metadata_env

    printf 'commit:  %s\n' "${AWMD_GIT_COMMIT:-<empty>}"
    printf 'branch:  %s\n' "${AWMD_GIT_BRANCH:-<empty>}"
    printf 'dirty:   %s\n' "${AWMD_GIT_DIRTY:-<empty>}"
    printf 'backend: %s\n' "${AWMD_BACKEND:-<empty>}"

    cd -- "${PROJECT_ROOT}"

    exec "${HALMD_CMD_ARRAY[@]}" "${BACKEND_ARGS[@]}"
}

main "$@"
