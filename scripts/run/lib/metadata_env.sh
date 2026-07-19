#!/usr/bin/env bash

set_metadata_env() {
    local repo="${PROJECT_ROOT:-.}"

    AWMD_GIT_COMMIT="$(
        git -C "${repo}" rev-parse HEAD 2>/dev/null || true
    )"

    AWMD_GIT_BRANCH="$(
        git -C "${repo}" branch --show-current 2>/dev/null || true
    )"

    if [[ -n "$(git -C "${repo}" status --porcelain 2>/dev/null)" ]]; then
        AWMD_GIT_DIRTY="true"
    else
        AWMD_GIT_DIRTY="false"
    fi

    AWMD_CREATED_AT="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

    AWMD_BACKEND="${AWMD_BACKEND:-native}"
    AWMD_DOCKER_IMAGE="${AWMD_DOCKER_IMAGE:-}"

    export AWMD_GIT_COMMIT
    export AWMD_GIT_BRANCH
    export AWMD_GIT_DIRTY
    export AWMD_CREATED_AT
    export AWMD_BACKEND
    export AWMD_DOCKER_IMAGE
}
