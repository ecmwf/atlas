#!/usr/bin/env bash
#
# (C) Copyright 2026- ECMWF.
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
#   --npromas "LIST"      Space-separated nproma values. Default depends on target:
#                         host: "8 16 32"
#                         device: "32 64 8192"
#   --nlevs "LIST"        Space-separated nlev values. Default: "1 32 137".
#   --nvars "LIST"        Space-separated nvar values (0 omits the variable dimension). Default: "0".
#   --nproma-other N      nproma for the "other" blocked field (b2b). Default: 64.
#   --npts N              Horizontal points. Default: 100000.
#   --iterations N        Timed iterations. Default: 50.
#   --warmup N            Warmup iterations. Default: 10.
#   --loop-orders "LIST"  Space-separated loop-orders to compare.
#                         Default: "nproma_innermost nproma_outermost".
#   --implementations "LIST" Space-separated implementations to compare.
#                         Default depends on target:
#                         host: "raw_pointers mdspan"
#                         device: "arrayview mdspan"
#   --on-device           Run on device (default).
#   --on-host             Run on host instead of device.
#   -v, --verbose         Print the launched benchmark command and stream its output.
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
PRECISIONS="float"
NLEVS="0 137"
NVARS="0 5 32"
NPTS=100000
ITERATIONS=50
WARMUP=10
ON_DEVICE=0
CSV=0
VERBOSE=0

HOST_DEFAULT_NPROMAS="8 16 32"
DEVICE_DEFAULT_NPROMAS="32 64 8192"
HOST_DEFAULT_NPROMA_OTHER=64
DEVICE_DEFAULT_NPROMA_OTHER=64
HOST_DEFAULT_LOOP_ORDERS="nproma_innermost nproma_outermost"
DEVICE_DEFAULT_LOOP_ORDERS="coalesced_read coalesced_write"
HOST_DEFAULT_IMPLEMENTATIONS="raw_pointers mdspan"
DEVICE_DEFAULT_IMPLEMENTATIONS="arrayview mdspan"

NPROMAS=""
NPROMA_OTHER=""
LOOP_ORDERS=""
IMPLEMENTATIONS=""
NPROMAS_SET=0
NPROMA_OTHER_SET=0
LOOP_ORDERS_SET=0
IMPLEMENTATIONS_SET=0

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
        --npromas)       NPROMAS="$2"; NPROMAS_SET=1; shift 2 ;;
        --nlevs)         NLEVS="$2"; shift 2 ;;
        --nvars)         NVARS="$2"; shift 2 ;;
        --nproma-other)  NPROMA_OTHER="$2"; NPROMA_OTHER_SET=1; shift 2 ;;
        --npts)          NPTS="$2"; shift 2 ;;
        --iterations)    ITERATIONS="$2"; shift 2 ;;
        --warmup)        WARMUP="$2"; shift 2 ;;
        --loop-orders)   LOOP_ORDERS="$2"; LOOP_ORDERS_SET=1; shift 2 ;;
        --implementations) IMPLEMENTATIONS="$2"; IMPLEMENTATIONS_SET=1; shift 2 ;;
        --on-device)     ON_DEVICE=1; shift ;;
        --on-host)       ON_DEVICE=0; shift ;;
        -v|--verbose)    VERBOSE=1; shift ;;
        --csv)           CSV=1; shift ;;
        -h|--help)       usage 0 ;;
        *) echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

if [[ "$ON_DEVICE" -eq 1 ]]; then
    if [[ "$NPROMAS_SET" -eq 0 ]]; then
        NPROMAS="$DEVICE_DEFAULT_NPROMAS"
    fi
    if [[ "$NPROMA_OTHER_SET" -eq 0 ]]; then
        NPROMA_OTHER="$DEVICE_DEFAULT_NPROMA_OTHER"
    fi
    if [[ "$LOOP_ORDERS_SET" -eq 0 ]]; then
        LOOP_ORDERS="$DEVICE_DEFAULT_LOOP_ORDERS"
    fi
    if [[ "$IMPLEMENTATIONS_SET" -eq 0 ]]; then
        IMPLEMENTATIONS="$DEVICE_DEFAULT_IMPLEMENTATIONS"
    fi
else
    if [[ "$NPROMAS_SET" -eq 0 ]]; then
        NPROMAS="$HOST_DEFAULT_NPROMAS"
    fi
    if [[ "$NPROMA_OTHER_SET" -eq 0 ]]; then
        NPROMA_OTHER="$HOST_DEFAULT_NPROMA_OTHER"
    fi
    if [[ "$LOOP_ORDERS_SET" -eq 0 ]]; then
        LOOP_ORDERS="$HOST_DEFAULT_LOOP_ORDERS"
    fi
    if [[ "$IMPLEMENTATIONS_SET" -eq 0 ]]; then
        IMPLEMENTATIONS="$HOST_DEFAULT_IMPLEMENTATIONS"
    fi
fi

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

# Colour support: bold escape sequences if stdout is a terminal that supports it.
if [[ -t 1 ]] && tput bold >/dev/null 2>&1; then
    BOLD="$(tput bold)"
    RESET="$(tput sgr0)"
else
    BOLD=""
    RESET=""
fi

print_command() {
    local arg
    for arg in "$@"; do
        printf '%q ' "$arg"
    done
    printf '\n'
}

# Run one benchmark configuration and echo lines "direction gelems" for each of the
# three relayout directions, parsed from the table output.
run_one() {
    local precision="$1" nproma="$2" nproma_other="$3" nlev="$4" nvar="$5" loop_order="$6" operation="$7" implementation="$8"
    local output output_file exit_code
    local CMD=(
        "$BIN"
        --npts="$NPTS" --nlev="$nlev" --nvar="$nvar"
        --nproma="$nproma" --nproma-other="$nproma_other"
        --iterations="$ITERATIONS" --warmup="$WARMUP"
        --precision="$precision" --loop-order="$loop_order"
        --operation="$operation" --implementation="$implementation" --format=table
        --blocked-nonblocked-use-memcpy
    )

    if [[ -n "$DEVICE_FLAG" ]]; then
        CMD+=("$DEVICE_FLAG")
    fi

    output_file="$(mktemp "${TMPDIR:-/tmp}/atlas-benchmark-relayout-output.XXXXXX")" || {
        echo "error: could not create temporary file for benchmark output" >&2
        return 1
    }

    # shellcheck disable=SC2086
    if [[ "$VERBOSE" -eq 1 ]]; then
        printf "=========================================================================================================================\n" >&2
        printf 'Running benchmark for precision=%s nproma=%s nproma-other=%s nlev=%s nvar=%s loop-order=%s operation=%s\n' \
            "$precision" "$nproma" "$nproma_other" "$nlev" "$nvar" "$loop_order" "$operation" >&2
        printf 'Command: ' >&2
        if [[ -n "$LAUNCH" ]]; then
            printf '%s ' "$LAUNCH" >&2
        fi
        print_command "${CMD[@]}" >&2

        if $LAUNCH "${CMD[@]}" 2>&1 | tee "$output_file" >&2
        then
            exit_code=0
        else
            exit_code=$?
        fi
    else
        if output="$($LAUNCH "${CMD[@]}" 2>&1)"
        then
            exit_code=0
            printf '%s\n' "$output" > "$output_file"
        else
            exit_code=$?
        fi
    fi

    output="$(cat "$output_file")"
    rm -f "$output_file"

    if [[ "$exit_code" -ne 0 ]]; then
        echo "error: benchmark failed for precision=$precision nproma=$nproma nproma-other=$nproma_other nlev=$nlev nvar=$nvar loop-order=$loop_order operation=$operation" >&2
        if [[ "$VERBOSE" -eq 0 ]]; then
            echo "$output" >&2
        fi
        return 1
    fi

    # Table rows: "<name> <min> <max> <avg> <stddev> <GB/s> <Gelem/s>".
    # Extract the Gelem/s (7th field) for each direction row.
    local d
    for d in $DIRECTIONS; do
        awk -v name="$d" '$1==name { print name, $7 }' <<<"$output"
    done
}

report_case() {
    local precision="$1" nproma="$2" nproma_other="$3" nlev="$4" nvar="$5" operation="$6" directions="$7"
    local gbs_file loop_order dir value best_lo best_impl best_val val key implementation

    gbs_file="$(mktemp "${TMPDIR:-/tmp}/atlas-benchmark-relayout.XXXXXX")" || {
        echo "error: could not create temporary file for benchmark results" >&2
        return 1
    }

    for implementation in $IMPLEMENTATIONS; do
        for loop_order in $LOOP_ORDERS; do
            while read -r dir value; do
                [[ -z "${dir:-}" ]] && continue
                printf '%s %s\n' "$dir|$loop_order|$implementation" "$value" >> "$gbs_file"
                if [[ "$CSV" -eq 1 ]]; then
                    echo "$DEVICE_LABEL,$precision,$NPTS,$nlev,$nvar,$nproma,$nproma_other,$implementation,$dir,$loop_order,$value"
                fi
            done < <(run_one "$precision" "$nproma" "$nproma_other" "$nlev" "$nvar" "$loop_order" "$operation" "$implementation")
        done
    done

    if [[ "$CSV" -eq 0 ]]; then
        local colw label ncols sepwidth sep n_impl n_lo header
        n_impl=0
        for implementation in $IMPLEMENTATIONS; do
            n_impl=$((n_impl + 1))
        done
        n_lo=0
        for loop_order in $LOOP_ORDERS; do
            n_lo=$((n_lo + 1))
        done

        colw=6
        ncols=0
        for implementation in $IMPLEMENTATIONS; do
            for loop_order in $LOOP_ORDERS; do
                if [[ "$n_impl" -eq 1 ]]; then
                    label="${loop_order}"
                elif [[ "$n_lo" -eq 1 ]]; then
                    label="${implementation}"
                else
                    label="${implementation}/${loop_order}"
                fi
                if [[ "${#label}" -gt "$colw" ]]; then
                    colw="${#label}"
                fi
                ncols=$((ncols + 1))
            done
        done
        # direction (2+24) + data/winner columns (ncols+1)*(1+colw) + speedup (1+8)
        sepwidth=$((26 + (ncols + 1) * (colw + 1) + 9))
        sep="$(printf '%*s' "$sepwidth" '' | tr ' ' '=')"

        header="target=$DEVICE_LABEL precision=$precision npts=$NPTS nlev=$nlev nvar=$nvar nproma=$nproma nproma_other=$nproma_other operation=$operation"
        if [[ "$n_impl" -eq 1 || "$n_lo" -eq 1 ]]; then
            header="target==$DEVICE_LABEL precision=$precision npts=$NPTS nlev=$nlev nvar=$nvar nproma=$nproma nproma_other=$nproma_other operation=$operation"
            if [[ "$n_impl" -eq 1 ]]; then
                header+=" implementation=$IMPLEMENTATIONS"
            fi
            if [[ "$n_lo" -eq 1 ]]; then
                header+=" loop_order=$LOOP_ORDERS"
            fi
        fi

        echo "$sep"
        echo "$header"
        printf '%*s\n' "$sepwidth" '' | tr ' ' '-'
        printf "  %-24s" "direction"
        for implementation in $IMPLEMENTATIONS; do
            for loop_order in $LOOP_ORDERS; do
                if [[ "$n_impl" -eq 1 ]]; then
                    label="${loop_order}"
                elif [[ "$n_lo" -eq 1 ]]; then
                    label="${implementation}"
                else
                    label="${implementation}/${loop_order}"
                fi
                printf " %*s" "$colw" "$label"
            done
        done
        printf " %*s %8s\n" "$colw" "winner" "speedup"
        for dir in $directions; do
            printf "  %-24s" "$dir"
            best_lo=""
            best_impl=""
            best_val=""
            worst_val=""
            for implementation in $IMPLEMENTATIONS; do
                for loop_order in $LOOP_ORDERS; do
                    key="$dir|$loop_order|$implementation"
                    val="$(awk -v key="$key" '$1==key { print $2; exit }' "$gbs_file")"
                    if [[ -n "$val" ]]; then
                        printf " %*s" "$colw" "$val"
                        if [[ -z "$best_val" ]] || awk -v a="$val" -v b="$best_val" 'BEGIN{exit !(a>b)}'; then
                            best_val="$val"
                            best_lo="$loop_order"
                            best_impl="$implementation"
                        fi
                        if [[ -z "$worst_val" ]] || awk -v a="$val" -v b="$worst_val" 'BEGIN{exit !(a<b)}'; then
                            worst_val="$val"
                        fi
                    else
                        printf " %*s" "$colw" "-"
                    fi
                done
            done
            local speedup pct speedup_plain winner winner_plain
            if [[ -n "$best_val" && -n "$worst_val" ]] && awk -v b="$worst_val" 'BEGIN{exit !(b>0)}'; then
                pct="$(awk -v best="$best_val" -v worst="$worst_val" 'BEGIN{printf "%.1f", (best-worst)/worst*100}')"
                speedup_plain="+${pct}%"
                speedup="$(printf "%8s" "$speedup_plain")"
                if [[ "$n_impl" -eq 1 ]]; then
                    winner_plain="${best_lo}"
                elif [[ "$n_lo" -eq 1 ]]; then
                    winner_plain="${best_impl}"
                else
                    winner_plain="${best_impl}/${best_lo}"
                fi
                winner="$(printf "%*s" "$colw" "$winner_plain")"
                if awk -v p="$pct" 'BEGIN{exit !(p>10)}'; then
                    speedup="${BOLD}${speedup}${RESET}"
                    winner="${BOLD}${winner}${RESET}"
                fi
            else
                speedup_plain="-"
                speedup="$(printf "%8s" "$speedup_plain")"
                winner="$(printf "%*s" "$colw" "n/a")"
            fi
            printf " %s %s\n" "$winner" "$speedup"
        done
        echo
    fi

    rm -f "$gbs_file"
}

# ---- sweep ------------------------------------------------------------------

if [[ "$CSV" -eq 1 ]]; then
    echo "target,precision,npts,nlev,nvar,nproma,nproma_other,implementation,direction,loop_order,gelems"
fi

for precision in $PRECISIONS; do
  for nproma in $NPROMAS; do
    for nlev in $NLEVS; do
      for nvar in $NVARS; do
        report_case "$precision" "$nproma" "$NPROMA_OTHER" "$nlev" "$nvar" "all" "$DIRECTIONS"

        if [[ "$ON_DEVICE" -eq 1 ]]; then
            report_case "$precision" "64" "8192" "$nlev" "$nvar" "b2b" "blocked_to_blocked"
        fi
     done
    done
  done
done
