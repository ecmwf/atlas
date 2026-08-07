#!/usr/bin/env python3

import argparse
from pathlib import Path


def read_power_spectrum(path):
    metadata = {}
    rows = []

    with path.open("r", encoding="utf-8") as file:
        for line in file:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                comment = stripped[1:].strip()
                if ":" in comment:
                    key, value = comment.split(":", 1)
                    metadata[key.strip().lower()] = value.strip()
                continue
            rows.append([float(value) for value in stripped.split()])

    if not rows:
        raise ValueError(f"No spectrum data found in {path}")

    columns = len(rows[0])
    if columns < 2:
        raise ValueError(f"Expected at least two columns in {path}")

    for row in rows:
        if len(row) != columns:
            raise ValueError(f"Inconsistent number of columns in {path}")

    wavenumbers = [int(row[0]) for row in rows]
    spectra_by_level = [[row[level_column] for row in rows] for level_column in range(1, columns)]
    return metadata, wavenumbers, spectra_by_level


def parse_args():
    parser = argparse.ArgumentParser(description="Plot an Atlas spectral power spectrum ASCII file.")
    parser.add_argument("spectrum", type=Path, help="Power-spectrum file written by atlas-icon-filter")
    parser.add_argument(
        "-l",
        "--level",
        type=int,
        default=1,
        help="1-based level index to plot when the spectrum contains multiple levels (default: 1)",
    )
    parser.add_argument("-o", "--output", type=Path, help="Write the plot to this image file instead of showing it")
    parser.add_argument("--title", help="Override the plot title")
    parser.add_argument("--logx", action="store_true", help="Use a logarithmic x-axis")
    parser.add_argument("--logy", action="store_true", help="Use a logarithmic y-axis")
    return parser.parse_args()


def main():
    args = parse_args()
    metadata, wavenumbers, spectra_by_level = read_power_spectrum(args.spectrum)

    nlev = len(spectra_by_level)
    if args.level < 1 or args.level > nlev:
        raise SystemExit(f"Requested level {args.level}, but {args.spectrum} contains {nlev} level(s)")

    import matplotlib.pyplot as plt

    power_spectrum = spectra_by_level[args.level - 1]
    field_name = metadata.get("field", args.spectrum.stem)
    title = args.title or f"Power spectrum: {field_name}, level {args.level}"

    fig, ax = plt.subplots()
    ax.plot(wavenumbers, power_spectrum, marker=".", linewidth=1.0)
    ax.set_xlabel("Wave number")
    ax.set_ylabel("Power spectrum")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)

    if args.logx:
        ax.set_xscale("log")
    if args.logy:
        ax.set_yscale("log")

    fig.tight_layout()
    if args.output:
        fig.savefig(args.output, dpi=160)
    else:
        plt.show()


if __name__ == "__main__":
    main()
