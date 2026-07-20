#!/usr/bin/env bash

require_h5ls() {
    if ! command -v h5ls >/dev/null 2>&1; then
        printf 'Error: h5ls is required to validate HDF5 output\n' >&2
        return 1
    fi
}

h5_is_readable() {
    local path="$1"
    [[ -f "${path}" ]] && h5ls "${path}" >/dev/null 2>&1
}

require_readable_h5() {
    local path="$1"
    if ! h5_is_readable "${path}"; then
        printf 'Error: missing or unreadable HDF5 file: %s\n' "${path}" >&2
        return 1
    fi
}
