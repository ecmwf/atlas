#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="$(cd -- "${ATLAS_BUILD_DIR:-$script_dir/../../../../../build-ninja}" && pwd)"
ninja_bin="$(sed -n 's/^CMAKE_MAKE_PROGRAM:FILEPATH=//p' "$build_dir/CMakeCache.txt")"
mpi_exec="$(sed -n 's/^MPIEXEC_EXECUTABLE:FILEPATH=//p' "$build_dir/CMakeCache.txt")"
mpi_numproc_flag="$(sed -n 's/^MPIEXEC_NUMPROC_FLAG:STRING=//p' "$build_dir/CMakeCache.txt")"
executable="$build_dir/atlas/src/sandbox/interpolation/atlas-global-unmatched-matrix"
interpolations="$build_dir/bin/atlas-interpolations"
source_grid="${SOURCE_GRID:-eORCA12_T}"
target_grid="${TARGET_GRID:-O160}"
interpolation="${INTERPOLATION:-finite-element}"
mpi_ranks="${MPI_RANKS:-2}"

if [[ -z "$ninja_bin" || ! -x "$ninja_bin" ]]; then
    echo "Ninja executable not found in $build_dir/CMakeCache.txt" >&2
    exit 127
fi

if [[ -z "$mpi_exec" || ! -x "$mpi_exec" || -z "$mpi_numproc_flag" ]]; then
    echo "MPI launcher not found in $build_dir/CMakeCache.txt" >&2
    exit 127
fi

if [[ ! "$mpi_ranks" =~ ^[1-9][0-9]*$ ]]; then
    echo "MPI_RANKS must be a positive integer" >&2
    exit 2
fi

"$ninja_bin" -C "$build_dir" atlas-interpolations atlas-global-unmatched-matrix

export ATLAS_ORCA_CACHING=1
export ATLAS_TRACE="${ATLAS_TRACE:-0}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"

run_dir="$(mktemp -d "$build_dir/unmatched-matrix.XXXXXX")"
echo "Output directory: $run_dir"
cd "$run_dir"
matrix=remap_matrix

"$interpolations" \
    --s.grid "$source_grid" \
    --t.grid "$target_grid" \
    --i.type "$interpolation" \
    --output-matrix \
    --matrix.name "$matrix" \
    --matrix.format eckit

if [[ ! -s "$matrix.eckit" ]]; then
    echo "Interpolation matrix was not created: $run_dir/$matrix.eckit" >&2
    exit 1
fi

exec "$mpi_exec" "$mpi_numproc_flag" "$mpi_ranks" "$executable" \
    --sgrid "$source_grid" \
    --tgrid "$target_grid" \
    --interpolation "$interpolation" \
    --matrix "$matrix" \
    --format eckit
