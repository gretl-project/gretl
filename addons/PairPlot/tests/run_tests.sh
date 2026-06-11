#!/usr/bin/env bash

set -u

readonly DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT="PairPlot"
readonly SCRIPTS=(
    "run_tests.inp"
    "pngtests.inp"
)

failures=0

cd -- "$DIR" || {
    printf 'Error: Cannot enter test directory: %s\n' "$DIR" >&2
    exit 1
}

printf 'Running tests for %s in %s\n\n' "$PROJECT" "$DIR"

if ! command -v gretlcli >/dev/null 2>&1; then
    printf 'Error: gretlcli was not found in PATH.\n' >&2
    exit 127
fi

for script in "${SCRIPTS[@]}"; do
    if [[ ! -f "$script" ]]; then
        printf 'Error: Test script not found: %s\n' "$script" >&2
        ((failures += 1))
        continue
    fi

    printf 'Running %s...\n' "$script"

    if gretlcli -b -e -q "$script"; then
        printf 'Success: %s passed.\n\n' "$script"
    else
        status=$?
        printf 'Failure: %s exited with status %d.\n\n' \
            "$script" "$status" >&2
        ((failures += 1))
    fi
done

if ((failures == 0)); then
    printf 'Success: All tests passed for %s.\n' "$PROJECT"
    exit 0
fi

printf 'Failure: %d test script(s) failed for %s.\n' \
    "$failures" "$PROJECT" >&2
exit 1