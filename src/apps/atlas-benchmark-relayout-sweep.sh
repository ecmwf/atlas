#!/usr/bin/env bash
#
# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
#
# atlas-benchmark-relayout-sweep.sh
#
# Sweep the atlas-benchmark-relayout loop-orders over a grid of precisions,
# nproma and nlev values, and report which loop-order achieves the highest
# memory bandwidth for each relayout direction (blocked_to_nonblocked,
# nonblocked_to_blocked, blocked_to_blocked).
#
# This is intended to characterize a given GPU: the best loop-order is device-,
# precision-, and shape-dependent (see the loop-order discussion in
# relayout_on_device.tcc), so sweeping and reading off the winners is more
# reliable than assuming a fixed default.
#
# Usage:
#   atlas-benchmark-relayout-sweep.sh [options]
#
# Options (all optional; sensible defaults are shown):
#   --bin PATH            Path to the atlas-benchmark-relayout binary.
#                         Default: $ATLAS_BENCHMARK_RELAYOUT, else ./bin/atlas-benchmark-relayout,
#                         else the first match found on $PATH.
#   --launch "CMD"        Launcher prefix used for every benchmark invocation, e.g.
#                         "srun -q dg --gpus-per-task=1 -n 1 -c 16". Default: none (run directly).
#                         Pass --launch "" to force a direct local run.
#   --precisions "LIST"   Space-separated value types. Default: "double float".
#   --npromas "LIST"      Space-separated nproma values. Default: "16 32 64".
#   --nlevs "LIST"        Space-separated nlev values. Default: "1 32 137".
#   --nvars "LIST"        Space-separated nvar values (0 omits the variable dimension). Default: "0".
#   --nproma-other N      nproma for the "other" blocked field (b2b). Default: 64.
#   --npts N              Horizontal points. Default: 100000.
#   --iterations N        Timed iterations. Default: 50.
#   --warmup N            Warmup iterations. Default: 10.
#   --loop-orders "LIST"  Space-separated loop-orders to compare.
#                         Default: "coalesced_write coalesced_read".
#   --on-device           Run on device (default).
#   --on-host             Run on host instead of device.
#   --csv                 Emit machine-readable CSV instead of the human-readable report.
#   -h, --help            Show this help and exit.
#
# Environment:
#   Source the build environment (e.g. env.sh) before running, so that the
#   benchmark binary and any GPU launcher are available.

set -euo pipefail

# ---- defaults ---------------------------------------------------------------
BIN="${ATLAS_BENCHMARK_RELAYOUT:-}"
LAUNCH=""
PRECISIONS="double float"
NPROMAS="16 32 64"
NLEVS="1 32 137"
NVARS="0"
NPROMA_OTHER=64
NPTS=100000
ITERATIONS=50
WARMUP=10
LOOP_ORDERS="coalesced_write coalesced_read"
ON_DEVICE=1
CSV=0

DIRECTIONS="blocked_to_nonblocked nonblocked_to_blocked blocked_to_blocked"

usage() {
    sed -n '2,64p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

# ---- argument parsing -------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bin)           BIN="$2"; shift 2 ;;
        --launch)        LAUNCH="$2"; shift 2 ;;
        --precisions)    PRECISIONS="$2"; shift 2 ;;
        --npromas)       NPROMAS="$2"; shift 2 ;;
        --nlevs)         NLEVS="$2"; shift 2 ;;
        --nvars)         NVARS="$2"; shift 2 ;;
        --nproma-other)  NPROMA_OTHER="$2"; shift 2 ;;
        --npts)          NPTS="$2"; shift 2 ;;
        --iterations)    ITERATIONS="$2"; shift 2 ;;
        --warmup)        WARMUP="$2"; shift 2 ;;
        --loop-orders)   LOOP_ORDERS="$2"; shift 2 ;;
        --on-device)     ON_DEVICE=1; shift ;;
        --on-host)       ON_DEVICE=0; shift ;;
        --csv)           CSV=1; shift ;;
        -h|--help)       usage 0 ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

# ---- locate the benchmark binary -------------------------------------------
if [[ -z "$BIN" ]]; then
    if [[ -x "./bin/atlas-benchmark-relayout" ]]; then
        BIN="./bin/atlas-benchmark-relayout"
    elif command -v atlas-benchmark-relayout >/dev/null 2>&1; then
        BIN="$(command -v atlas-benchmark-relayout)"
    fi
fi
if [[ -z "$BIN" || ! -x "$BIN" ]]; then
    echo "error: could not find the atlas-benchmark-relayout binary." >&2
    echo "       pass --bin PATH or set ATLAS_BENCHMARK_RELAYOUT." >&2
    exit 1
fi

DEVICE_FLAG=""
DEVICE_LABEL="host"
if [[ "$ON_DEVICE" -eq 1 ]]; then
    DEVICE_FLAG="--on-device"
    DEVICE_LABEL="device"
fi

# ---- helpers ----------------------------------------------------------------

# Run one benchmark configuration and echo lines "direction gbs" for each of the
# three relayout directions, parsed from the table output.
run_one() {
    local precision="$1" nproma="$2" nlev="$3" nvar="$4" loop_order="$5"
    local output
    # shellcheck disable=SC2086
    output="$($LAUNCH "$BIN" \
        --npts="$NPTS" --nlev="$nlev" --nvar="$nvar" \
        --nproma="$nproma" --nproma-other="$NPROMA_OTHER" \
        --iterations="$ITERATIONS" --warmup="$WARMUP" \
        --precision="$precision" --loop-order="$loop_order" \
        --format=table $DEVICE_FLAG 2>&1)" || {
        echo "error: benchmark failed for precision=$precision nproma=$nproma nlev=$nlev nvar=$nvar loop-order=$loop_order" >&2
        echo "$output" >&2
        return 1
    }
    # Table rows: "<name> <min> <max> <avg> <stddev> <GB/s> <Gelem/s>".
    # Extract the GB/s (6th field) for each direction row.
    local d
    for d in $DIRECTIONS; do
        awk -v name="$d" '$1==name { print name, $6 }' <<<"$output"
    done
}

# ---- sweep ------------------------------------------------------------------

if [[ "$CSV" -eq 1 ]]; then
    echo "target,precision,npts,nlev,nvar,nproma,nproma_other,direction,loop_order,gbs"
fi

for precision in $PRECISIONS; do
  for nproma in $NPROMAS; do
    for nlev in $NLEVS; do
      for nvar in $NVARS; do

        # Collect GB/s per direction per loop-order into associative arrays.
        declare -A gbs=()
        for loop_order in $LOOP_ORDERS; do
            while read -r dir value; do
                [[ -z "${dir:-}" ]] && continue
                gbs["$dir|$loop_order"]="$value"
                if [[ "$CSV" -eq 1 ]]; then
                    echo "$DEVICE_LABEL,$precision,$NPTS,$nlev,$nvar,$nproma,$NPROMA_OTHER,$dir,$loop_order,$value"
                fi
            done < <(run_one "$precision" "$nproma" "$nlev" "$nvar" "$loop_order")
        done

        if [[ "$CSV" -eq 0 ]]; then
            echo "===================================================================================="
            echo "target=$DEVICE_LABEL precision=$precision npts=$NPTS nlev=$nlev nvar=$nvar nproma=$nproma nproma_other=$NPROMA_OTHER"
            echo "------------------------------------------------------------------------------------"
            printf "  %-24s" "direction"
            for loop_order in $LOOP_ORDERS; do
                printf " %18s" "$loop_order"
            done
            printf " %18s\n" "winner"
            for dir in $DIRECTIONS; do
                printf "  %-24s" "$dir"
                best_lo=""
                best_val=""
                for loop_order in $LOOP_ORDERS; do
                    val="${gbs["$dir|$loop_order"]:-}"
                    if [[ -n "$val" ]]; then
                        printf " %18s" "$val"
                        if [[ -z "$best_val" ]] || awk -v a="$val" -v b="$best_val" 'BEGIN{exit !(a>b)}'; then
                            best_val="$val"
                            best_lo="$loop_order"
                        fi
                    else
                        printf " %18s" "-"
                    fi
                done
                printf " %18s\n" "${best_lo:-n/a}"
            done
            echo
        fi

        unset gbs
     done
    done
  done
done
