#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd -- "${SCRIPT_DIR}/../../.." && pwd)}"

main() {
    local image
    local name
    local user_name
    local -a docker_args=(run --rm)

    if [[ -z "${AWMD_DOCKER_IMAGE:-}" ]]; then
        printf 'Error: AWMD_DOCKER_IMAGE is not set\n' >&2
        return 2
    fi
    image="${AWMD_DOCKER_IMAGE}"

    if [[ ! -d "${PROJECT_ROOT}" ]]; then
        printf 'Error: project root not found: %s\n' "${PROJECT_ROOT}" >&2
        return 1
    fi

    if [[ -t 0 && -t 1 ]]; then
        docker_args+=(-it)
    fi

    if [[ "${AWMD_DOCKER_GPUS:-all}" != "none" ]]; then
        docker_args+=(--gpus "${AWMD_DOCKER_GPUS:-all}")
    fi

    if [[ "${RUN_AS_ROOT:-0}" != "1" ]]; then
        user_name="$(id -un)"
        docker_args+=(
            --user "$(id -u):$(id -g)"
            --mount "type=bind,src=/etc/passwd,dst=/etc/passwd,readonly"
            --mount "type=bind,src=/etc/group,dst=/etc/group,readonly"
            --env "USER=${user_name}"
            --env "LOGNAME=${user_name}"
        )
    fi

    while IFS= read -r name; do
        [[ "${name}" == AWMD_* ]] || continue
        docker_args+=(--env "${name}")
    done < <(compgen -e)

    docker_args+=(
        --mount "type=bind,src=${PROJECT_ROOT},dst=/work"
        --workdir /work
        --env HOME=/work
        "${image}"
    )
    docker_args+=("$@")

    exec docker "${docker_args[@]}"
}

main "$@"
