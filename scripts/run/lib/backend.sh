#!/usr/bin/env bash

parse_backend_args() {
    local backend="native"
    local custom_image=""

    BACKEND_ARGS=()

    while (( $# > 0 )); do
        case "$1" in
            --)
                shift
                BACKEND_ARGS+=("$@")
                break
                ;;

            --native)
                backend="native"
                shift
                ;;

            --docker)
                backend="docker"
                shift
                ;;

            --docker-position-shifting)
                backend="docker-position-shifting"
                shift
                ;;

            --docker-image)
                if (( $# < 2 )); then
                    printf 'Error: --docker-image requires an image name\n' >&2
                    return 2
                fi

                backend="docker-custom"
                custom_image="$2"
                shift 2
                ;;

            --docker-image=*)
                backend="docker-custom"
                custom_image="${1#*=}"
                shift
                ;;

            *)
                BACKEND_ARGS+=("$1")
                shift
                ;;
        esac
    done

    case "$backend" in
        native)
            AWMD_BACKEND="native"
            AWMD_DOCKER_IMAGE=""
            HALMD_CMD_ARRAY=(halmd)
            ;;

        docker)
            AWMD_BACKEND="docker"
            AWMD_DOCKER_IMAGE="halmd"
            HALMD_CMD_ARRAY=(
                "${PROJECT_ROOT}/scripts/run/docker/halmd.sh"
                halmd
            )
            ;;

        docker-position-shifting)
            AWMD_BACKEND="docker"
            AWMD_DOCKER_IMAGE="halmd:position-shifting"
            HALMD_CMD_ARRAY=(
                "${PROJECT_ROOT}/scripts/run/docker/halmd.sh"
                halmd
            )
            ;;

        docker-custom)
            if [[ -z "${custom_image}" ]]; then
                printf 'Error: custom Docker image cannot be empty\n' >&2
                return 2
            fi

            AWMD_BACKEND="docker"
            AWMD_DOCKER_IMAGE="${custom_image}"
            HALMD_CMD_ARRAY=(
                "${PROJECT_ROOT}/scripts/run/docker/halmd.sh"
                halmd
            )
            ;;
    esac

    export AWMD_BACKEND
    export AWMD_DOCKER_IMAGE
}
