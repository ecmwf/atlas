#!/usr/bin/env python3
# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Read and visualise interpolation stencils extracted by

    atlas-interpolations --output-stencil ...

The ASCII file written by ``--output-stencil`` contains one or more stencils.
Each stencil holds:
  * one target point (record ``T``),
  * the source points contributing to it, with weights (records ``S``),
  * extra source points around the stencil, for context (records ``X``),
  * extra target points around the target, for context (records ``Y``),
  * optionally the source-mesh cells touching the stencil (records ``C``),
  * optionally the target-mesh cells around the target (records ``D``).

This module provides a small parser (:func:`read_stencils`) returning plain
Python/NumPy data, and a simple matplotlib viewer (:func:`plot_stencil`).

Usage:
    ./plot_stencil.py stencil.txt                     # all stencils as subfigures in one plot
    ./plot_stencil.py stencil.txt --index 0           # only the stencil with index 0
    ./plot_stencil.py stencil.txt --save out.png      # save the figure to out.png
    ./plot_stencil.py stencil.txt --stereographic     # stereographic projection around each target
    ./plot_stencil.py stencil.txt --center=10,45      # stereographic projection around (lon,lat)=(10,45)
    ./plot_stencil.py stencil.txt --weight-color      # colour/size stencil by weight; mask on the edge

Each subfigure's title shows the target's 0-based global index in the target grid and its (lon, lat).
The parser has no dependencies; plotting requires matplotlib.
"""

import argparse
import math
import sys
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple


@dataclass
class Point:
    gidx: int
    lon: float
    lat: float
    mask: int
    weight: Optional[float] = None  # set for stencil (S) points only


@dataclass
class Stencil:
    index: int
    target: Optional[Point] = None
    stencil: List[Point] = field(default_factory=list)  # source points with weights (S)
    context: List[Point] = field(default_factory=list)  # context source points (X)
    target_context: List[Point] = field(default_factory=list)  # context target points (Y)
    cells: List[List[int]] = field(default_factory=list)  # source-mesh cells (source gidx lists, C)
    target_cells: List[List[int]] = field(default_factory=list)  # target-mesh cells (target gidx lists, D)

@dataclass
class StencilFile:
    meta: dict = field(default_factory=dict)
    stencils: List[Stencil] = field(default_factory=list)


def read_stencils(path: str) -> StencilFile:
    """Parse a stencil file written by ``atlas-interpolations --output-stencil``."""
    result = StencilFile()
    current: Optional[Stencil] = None

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                # metadata lines look like '# key=value'
                body = line[1:].strip()
                if "=" in body:
                    key, _, value = body.partition("=")
                    result.meta[key.strip()] = value.strip()
                continue

            tok = line.split()
            kind = tok[0]

            if kind == "BEGIN":
                current = Stencil(index=int(tok[1]))
            elif kind == "END":
                if current is not None:
                    result.stencils.append(current)
                current = None
            elif kind == "T":
                gidx, lon, lat, mask = int(tok[1]), float(tok[2]), float(tok[3]), int(tok[4])
                current.target = Point(gidx, lon, lat, mask)
            elif kind == "S":
                gidx, lon, lat = int(tok[1]), float(tok[2]), float(tok[3])
                weight, mask = float(tok[4]), int(tok[5])
                current.stencil.append(Point(gidx, lon, lat, mask, weight))
            elif kind == "X":
                gidx, lon, lat, mask = int(tok[1]), float(tok[2]), float(tok[3]), int(tok[4])
                current.context.append(Point(gidx, lon, lat, mask))
            elif kind == "Y":
                gidx, lon, lat, mask = int(tok[1]), float(tok[2]), float(tok[3]), int(tok[4])
                current.target_context.append(Point(gidx, lon, lat, mask))
            elif kind == "C":
                n = int(tok[1])
                current.cells.append([int(g) for g in tok[2:2 + n]])
            elif kind == "D":
                n = int(tok[1])
                current.target_cells.append([int(g) for g in tok[2:2 + n]])
            else:
                raise ValueError(f"Unrecognised record '{kind}' in {path}: {line!r}")

    return result


# A projector maps a geographic (lon, lat) in degrees to plane coordinates (x, y).
Projector = Callable[[float, float], Tuple[float, float]]

_RAD2DEG = 180.0 / math.pi


def lonlat_projector(ref_lon: float) -> Projector:
    """Plate-carree-like projection: lon/lat as-is, with lon unwrapped to within 180 deg of ref_lon."""
    def project(lon: float, lat: float) -> Tuple[float, float]:
        while lon - ref_lon > 180.0:
            lon -= 360.0
        while lon - ref_lon < -180.0:
            lon += 360.0
        return lon, lat
    return project


def stereographic_projector(lon0: float, lat0: float) -> Projector:
    """Azimuthal stereographic projection centred on (lon0, lat0).

    Handles the poles and the date line naturally and is conformal (locally shape-preserving),
    which makes it well suited to inspecting a stencil around a point. Plane units are scaled so
    that, near the centre, one unit is approximately one degree of great-circle arc.
    """
    lon0r = math.radians(lon0)
    lat0r = math.radians(lat0)
    sin_lat0 = math.sin(lat0r)
    cos_lat0 = math.cos(lat0r)

    def project(lon: float, lat: float) -> Tuple[float, float]:
        lonr = math.radians(lon)
        latr = math.radians(lat)
        dlon = lonr - lon0r
        cos_c = sin_lat0 * math.sin(latr) + cos_lat0 * math.cos(latr) * math.cos(dlon)
        denom = max(1.0 + cos_c, 1e-12)  # guard against the antipode
        k = 2.0 / denom
        x = k * math.cos(latr) * math.sin(dlon)
        y = k * (cos_lat0 * math.sin(latr) - sin_lat0 * math.cos(latr) * math.cos(dlon))
        return x * _RAD2DEG, y * _RAD2DEG
    return project


# --- Weight-colour mode helpers -----------------------------------------------------------------
# The mask normally uses colour (green = valid, red = masked). To colour points by weight instead
# without a clash, weight-colour mode uses a red-free sequential colormap (viridis) for the weight
# and moves the mask onto the marker edge: a neutral ring for valid, a bold red ring for masked.
WEIGHT_CMAP = "custom_greens"
# WEIGHT_CMAP = "jet"
# WEIGHT_CMAP = "Greens"

def _weights_cmap():
    if WEIGHT_CMAP != "custom_greens":
        import matplotlib.pyplot as plt
        return plt.get_cmap(WEIGHT_CMAP)
    else:
        from matplotlib.colors import LinearSegmentedColormap
        return LinearSegmentedColormap.from_list("custom_greens", ["#E8F5E9", "#81C784", "#1B5E20"])

def _weight_size(weight_size_dynamic: bool, frac: Optional[float]) -> float:
    """Marker area for a weight fraction in [0, 1] (redundant size encoding, like the reference)."""
    base_size = 40.0
    if weight_size_dynamic:
        f = 0.0 if frac is None else max(0.0, min(1.0, float(frac)))
        return base_size * (1. + 8.0*f)
    else:
        return 4.0*base_size

def _mask_edges(points, valid_edge: str):
    """Per-point edge colours / widths encoding the mask (red ring = masked)."""
    #edges = ["red" if not p.mask else valid_edge for p in points]
    #lws = [2.4 if not p.mask else 1.0 for p in points]
    edges = ["red" if not p.mask else valid_edge for p in points]
    lws = [1.0 if not p.mask else 1.0 for p in points]
    return edges, lws


def _boundary_lines(coords, mask):
    # =====================================================================
    # 2. COMPUTE BOUNDARY BLINDLY USING VORONOI
    # =====================================================================
    import numpy as np
    from scipy.spatial import Voronoi

    pts = np.array([coords[g] for g in coords])
    active_mask = [mask[g] for g in mask]

    # Calculate the physical bounds of your true node domain
    x_min, x_max = pts[:, 0].min(), pts[:, 0].max()
    y_min, y_max = pts[:, 1].min(), pts[:, 1].max()
    
    # Add a slight padding buffer (e.g., 5% of domain size) 
    # to avoid clipping valid boundary vertices right on the outer crust
    x_pad = (x_max - x_min) * 0.05
    y_pad = (y_max - y_min) * 0.05

    vor = Voronoi(pts)
    boundary_lines = []

    for ridge_idx, point_indices in enumerate(vor.ridge_points):
        p1_idx, p2_idx = point_indices
        
        # Extract only if the ridge sits between opposite flags (XOR logic) [1]
        if active_mask[p1_idx] != active_mask[p2_idx]:
            v_indices = vor.ridge_vertices[ridge_idx]
            
            if -1 not in v_indices:
                # Finite segment inside the point envelope
                v1, v2 = vor.vertices[v_indices]

                v1_inside = (x_min - x_pad <= v1[0] <= x_max + x_pad) and (y_min - y_pad <= v1[1] <= y_max + y_pad)
                v2_inside = (x_min - x_pad <= v2[0] <= x_max + x_pad) and (y_min - y_pad <= v2[1] <= y_max + y_pad)
                
                if v1_inside and v2_inside:
                    boundary_lines.append((v1, v2))
                else:
                    # Drop the segment if it shoots out of bounds
                    continue   
            else:
                continue


    return boundary_lines


def plot_stencil(sten: Stencil, meta: dict, project: Projector, ax=None,
                 xlabel: str = "lon [deg]", ylabel: str = "lat [deg]", show_legend: bool = True,
                 weight_color: bool = False, weight_size_dynamic: bool = False, wnorm=None, wcmap=None,
                 context: bool = False, cells: bool = True,
                 target_context: bool = False, target_cells: bool = False):
    """Plot a single stencil (source/target points, weights, masks, mesh cells) into `ax`.

    `project` maps geographic (lon, lat) to plane coordinates (see lonlat_projector /
    stereographic_projector). With `weight_color`, stencil points are coloured/sized by their
    weight (colormap `wcmap`, normalisation `wnorm`) and the mask is shown by the marker edge.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection, LineCollection

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    # Separate gidx -> (x, y) lookups for the source and target grids (their global
    # indices live in independent index spaces), so cells can be drawn from their gidx lists.
    src_coords = {p.gidx: project(p.lon, p.lat) for p in sten.context + sten.stencil}
    tgt_coords = {p.gidx: project(p.lon, p.lat) for p in sten.target_context + ([sten.target] if sten.target else [])}

    src_mask = {p.gidx: p.mask for p in sten.context + sten.stencil}
    tgt_mask = {p.gidx: p.mask for p in sten.target_context + ([sten.target] if sten.target else [])}

    combined_coords = {**src_coords, **tgt_coords}
    combined_mask = {**src_mask, **tgt_mask}
    src_mask_boundary_lines = _boundary_lines(combined_coords, combined_mask)

    def draw_cells(cells, coords, **kw):
        patches = []
        for cell in cells:
            pts = [coords[g] for g in cell if g in coords]
            if len(pts) == len(cell):
                patches.append(Polygon(pts, closed=True))
        if patches:
            ax.add_collection(PatchCollection(patches, **kw))

    # Mesh cells (drawn first, underneath the points): source solid grey, target dashed blue
    if target_cells:
        draw_cells(sten.target_cells, tgt_coords, facecolor="none", edgecolor="steelblue",
                   linewidths=0.6, linestyles="--", zorder=0)
    if cells:
        draw_cells(sten.cells, src_coords, facecolor="none", edgecolor="0.7", linewidths=0.6, zorder=1)

    def draw_points(points, marker, base_size, valid_edge, valid_lw, zorder):
        if not points:
            return
        xs = [project(p.lon, p.lat)[0] for p in points]
        ys = [project(p.lon, p.lat)[1] for p in points]
        if not weight_color:
            if points[0].weight is not None:  # stencil points carry a weight
                w = [p.weight for p in points]
                sizes = [_weight_size(weight_size_dynamic, wnorm(val)) for val in w]
            else:
                sizes = base_size
            # Default: colour encodes the mask (green = valid, red = masked)
            ax.scatter(xs, ys, c=[p.mask for p in points], cmap="RdYlGn", vmin=0, vmax=1,
                    s=sizes, marker=marker, edgecolors=valid_edge, linewidths=valid_lw,
                    zorder=zorder)
            return
        # Weight-colour mode: mask on the marker edge (red ring = masked)
        edges, lws = _mask_edges(points, valid_edge)
        if points[0].weight is not None:  # stencil points carry a weight
            w = [p.weight for p in points]
            sizes = [_weight_size(weight_size_dynamic, wnorm(val)) for val in w]
            ax.scatter(xs, ys, c=w, cmap=wcmap, norm=wnorm, s=sizes, marker=marker,
                       edgecolors=edges, linewidths=lws, zorder=zorder)
        else:  # context points have no weight -> neutral grey fill
            ax.scatter(xs, ys, c="0.8", s=base_size, marker=marker,
                       edgecolors=edges, linewidths=lws, zorder=zorder)

    # Context target points (diamonds), context source points, and weighted stencil points
    if target_context:
        draw_points(sten.target_context, "D", 55, "steelblue", 1.0, 2)
    if context:
        draw_points(sten.context, "o", 40, "0.4", 0.5, 2)
    draw_points(sten.stencil, "o", 160, "black", 1.2, 3)
    for p in sten.stencil:
        ax.annotate(
            f"{p.weight:.3g}", project(p.lon, p.lat),
            textcoords="offset points", xytext=(6, 6), fontsize=8, zorder=5,
        )

    # Target point
    if sten.target:
        tx, ty = project(sten.target.lon, sten.target.lat)
        ax.plot(tx, ty, marker="*", markersize=20, color="royalblue",
                markeredgecolor="black", zorder=4, label="target")

    if src_mask_boundary_lines:
        lc = LineCollection(src_mask_boundary_lines, colors='crimson', linewidths=1.0, zorder=4, label='Extracted Boundary Path')
        ax.add_collection(lc)

    # Title: 0-based global index in the target grid, and its coordinate
    if sten.target:
        gidx0 = sten.target.gidx - 1  # file stores 1-based; report 0-based as requested
        loc = f"target #{gidx0}  ({sten.target.lon:.4g}, {sten.target.lat:.4g})"
    else:
        loc = "target #?"
    ax.set_title(
        f"{loc}\n{meta.get('source_grid', '?')} -> {meta.get('target_grid', '?')}  "
        f"[{meta.get('interpolation', '?')}]",
        fontsize=10,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, linestyle=":", alpha=0.4)

    if show_legend:
        label = "stencil source (colour/size = weight)" if weight_color else "stencil source (weighted)"
        ax.legend(handles=_legend_handles(weight_color, label), fontsize=7, loc="best", framealpha=0.9)
    return ax


def _legend_handles(weight_color: bool, stencil_label: str = "stencil source (weighted)"):
    """Common legend proxy artists for both 2D and 3D, depending on the colour mode."""
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="*", color="w", markerfacecolor="royalblue",
               markeredgecolor="black", markersize=14, label="target"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="0.6",
               markeredgecolor="black", markersize=11, label=stencil_label),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="0.8",
               markeredgecolor="0.4", markersize=7, label="context source"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="0.8",
               markeredgecolor="steelblue", markersize=8, label="context target"),
    ]
    if weight_color:
        handles += [
            Line2D([0], [0], marker="o", color="w", markerfacecolor="0.8", markeredgecolor="0.3",
                   markersize=10, label="mask = 1 (valid)"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="0.8", markeredgecolor="red",
                   markeredgewidth=1.0, markersize=10, label="mask = 0 (masked, red ring)"),
        ]
    else:
        handles += [
            Line2D([0], [0], marker="s", color="w", markerfacecolor="tab:green",
                   markersize=9, label="mask = 1 (valid)"),
            Line2D([0], [0], marker="s", color="w", markerfacecolor="tab:red",
                   markersize=9, label="mask = 0 (masked)"),
        ]
    return handles


def _make_projector(sten: Stencil, stereographic: bool, center: Optional[Tuple[float, float]]):
    """Choose a projector and axis labels for one stencil."""
    if stereographic or center is not None:
        lon0, lat0 = center if center is not None else (
            (sten.target.lon, sten.target.lat) if sten.target else (0.0, 0.0))
        label = f"[deg from ({lon0:.4g}, {lat0:.4g})]"
        return stereographic_projector(lon0, lat0), f"x {label}", f"y {label}"
    ref_lon = sten.target.lon if sten.target else 0.0
    return lonlat_projector(ref_lon), "lon [deg]", "lat [deg]"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file", help="stencil file written by atlas-interpolations --output-stencil")
    parser.add_argument("--index", type=int, default=None, help="only plot the stencil with this index")
    parser.add_argument("--save", default=None, help="save the figure to this PNG path instead of showing")
    parser.add_argument("--stereographic", "--stereo", action="store_true", dest="stereographic",
                        help="use a stereographic projection centred on each stencil's target point")
    parser.add_argument("--center", "--centre", default=None, dest="center",
                        help="projection centre 'lon,lat' shared by all subplots (implies --stereographic)")
    parser.add_argument("--weight-color", "--weight-colour", action="store_true", dest="weight_color",
                        help="colour (and size) stencil points by their weight, with a colorbar; the "
                             "mask is then shown by the marker edge (red ring = masked)")
    parser.add_argument("--weight-size-dynamic", action="store_true", dest="weight_size_dynamic",
                        help="dynamically adjust the size of stencil points based on their weight, in addition to colouring them; the "
                             "mask is then shown by the marker edge (red ring = masked)")
    parser.add_argument("--context", action="store_true", help="show the context source points for each stencil")
    parser.add_argument("--cells", action="store_true", help="show the source cells for each stencil including context")
    parser.add_argument("--target-context", action="store_true", help="show the context points for the target of each stencil")
    parser.add_argument("--target-cells", action="store_true", help="show the target cells in the target context context")
    args = parser.parse_args(argv)

    center = None
    if args.center is not None:
        parts = args.center.split(",")
        if len(parts) != 2:
            print("--center expects 'lon,lat'", file=sys.stderr)
            return 1
        center = (float(parts[0]), float(parts[1]))

    data = read_stencils(args.file)
    print(f"# {args.file}")
    for k, v in data.meta.items():
        print(f"#   {k} = {v}")
    print(f"# parsed {len(data.stencils)} stencil(s)")

    selected = [s for s in data.stencils if args.index is None or s.index == args.index]
    if not selected:
        print(f"No stencil with index {args.index}", file=sys.stderr)
        return 1

    import matplotlib.pyplot as plt

    # Weight colour scale, shared across all subplots for a consistent colorbar.
    wnorm = wcmap = None
    if args.weight_color or args.weight_size_dynamic:
        from matplotlib.colors import Normalize, LinearSegmentedColormap
        all_w = [p.weight for s in selected for p in s.stencil if p.weight is not None]
        wmax = max(all_w) if all_w else 1.0
        wmin = min(all_w) if all_w else 0.0
        wnorm = Normalize(vmin=min(0.0, wmin), vmax=max(wmax, 1e-12))
        wcmap = _weights_cmap()

    # Arrange all selected stencils as subfigures in a single figure (near-square grid).
    n = len(selected)
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    subplot_kw = {}
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, 6.0 * nrows),
                             squeeze=False, subplot_kw=subplot_kw)
    axes_flat = [ax for row in axes for ax in row]

    plot = plot_stencil
    for i, (ax, sten) in enumerate(zip(axes_flat, selected)):
        project, xlabel, ylabel = _make_projector(sten, args.stereographic, center)
        plot(sten, data.meta, project, ax=ax, xlabel=xlabel, ylabel=ylabel,
             show_legend=(i == 0),  # one legend is enough for the whole figure
             weight_color=args.weight_color,
             weight_size_dynamic=args.weight_size_dynamic, wnorm=wnorm, wcmap=wcmap,
             context=args.context, cells=args.cells,
             target_context=args.target_context, target_cells=args.target_cells)
    for ax in axes_flat[n:]:
        ax.axis("off")  # hide unused cells in the grid

    fig.tight_layout()
    if args.weight_color:
        from matplotlib.cm import ScalarMappable
        sm = ScalarMappable(norm=wnorm, cmap=wcmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes_flat[:n], orientation="horizontal",
                            shrink=0.6, pad=0.08, aspect=40)
        cbar.set_label("interpolation weight")
    if args.save:
        out = args.save if args.save.lower().endswith(".png") else args.save + ".png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
