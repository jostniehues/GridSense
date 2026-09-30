import json
import sys

import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Patch, Rectangle
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["cmr10"],
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,
    "font.size": 8,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "figure.titlesize": 10,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 9,
    "axes.labelpad": 2,
    "figure.constrained_layout.w_pad": 0.005,
    "figure.constrained_layout.h_pad": 0.005,
})


def pdg_title(pdg):
    if pdg is None:
        return None
    pdg_dic = {
    0: "Gluon" ,   
    1: "Down-Quark",
    2: "Up-Quark",
    3: "Strange-Quark",
    4: "Charm-Quark",
    5: "Bottom-Quark",
    6: "Top-Quark",
    21: "Gluon",

    -1: "Anti-Down-Quark",
    -2: "Anti-Up-Quark",
    -3: "Anti-Strange-Quark",
    -4: "Anti-Charm-Quark",
    -5: "Anti-Bottom-Quark",
    -6: "Anti-Top-Quark",
}

    return f"{pdg_dic[pdg]}"

def as_title(order):
    if order is None:
        return f"All Orders"

    labels = {0: "LO", 1: "NLO", 2:"NNLO"}
    return labels.get(order, rf"$\mathrm{{N}}^{{{order}}}\mathrm{{LO}}$")
    



def support_boundary_segments(x_edges, q2_edges, selected):
    """Return the exposed cell edges of a boolean support mask."""
    segments = []
    nq, nx = selected.shape

    for iq in range(nq):
        for ix in range(nx):
            if not selected[iq, ix]:
                continue

            x_left = x_edges[ix]
            x_right = x_edges[ix + 1]
            q_bottom = q2_edges[iq]
            q_top = q2_edges[iq + 1]

            if ix == 0 or not selected[iq, ix - 1]:
                segments.append(((x_left, q_bottom), (x_left, q_top)))
            if ix == nx - 1 or not selected[iq, ix + 1]:
                segments.append(((x_right, q_bottom), (x_right, q_top)))
            if iq == 0 or not selected[iq - 1, ix]:
                segments.append(((x_left, q_bottom), (x_right, q_bottom)))
            if iq == nq - 1 or not selected[iq + 1, ix]:
                segments.append(((x_left, q_top), (x_right, q_top)))

    return segments


def single_color_map(color, name):
    """Build a light-to-full gradient from one process color."""
    rgb = np.asarray(to_rgb(color))
    #light_rgb = 1.0 - 0.50 * (1.0 - rgb)
    light_rgb = 1.0 - 0.15 * (1.0 - rgb)
    return LinearSegmentedColormap.from_list(
        f"{name}_gradient",
        [light_rgb, rgb],
    )


def draw_overlap(ax, grouped_support):
    """Hatch cells selected by at least two displayed support maps."""
    supports = list(grouped_support.values())
    if len(supports) < 2:
        return

    reference = supports[0]
    x_edges = reference["x_edges"]
    q2_edges = reference["q2_edges"]
    overlap_count = np.zeros_like(reference["selected"], dtype=np.uint16)

    for support in supports:
        if not (
            np.array_equal(x_edges, support["x_edges"])
            and np.array_equal(q2_edges, support["q2_edges"])
        ):
            return
        overlap_count += support["selected"]

    overlap = overlap_count >= 2
    for iq, ix in np.argwhere(overlap):
        ax.add_patch(
            Rectangle(
                (x_edges[ix], q2_edges[iq]),
                x_edges[ix + 1] - x_edges[ix],
                q2_edges[iq + 1] - q2_edges[iq],
                facecolor=(0.2, 0.2, 0.2, 0.10),
                edgecolor=(0.12, 0.12, 0.12, 0.30),
                linewidth=0.0,
                hatch="///",
                zorder=2.5,
            )
        )


def main() -> int:
    if len(sys.argv) != 6:
        print(
            "usage: python/plot_auxgrid.py "
            "<input.json> <output.pdf> <plot-style> <mode> <subtitle>",
            file=sys.stderr,
        )
        return 1

    input_json = sys.argv[1]
    output_pdf = sys.argv[2]
    plot_style = sys.argv[3]
    mode = sys.argv[4]
    subtitle = sys.argv[5]

    if plot_style not in {"boundary", "heatmap"}:
        print(f"unknown plot style: {plot_style}", file=sys.stderr)
        return 1

    if mode not in {"supp", "sense"}:
        print(f"unknown mode: {mode}", file=sys.stderr)
        return 1

    with open(input_json, encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, list) or len(data) == 0:
        print("input JSON does not contain any auxiliary grids", file=sys.stderr)
        return 1

    colors = [
        "#C87524",
        "#079B80",
        "#B65379",
        "#728F48",
        "#604A96",
        "#A53D50",
        "#916333",
        "#386F94",
    ]

    line_styles = [
        "--",
        "-.",
        (0, (5, 2, 1, 2)),
        (0, (3, 1)),
        (0, (1, 1)),
    ]

    pdg_groups = {}
    grid_colors = {}
    grid_linestyles = {}
    legend_handles = {}

    for entry in data:
        if not isinstance(entry, dict):
            print("input JSON entry has unexpected structure", file=sys.stderr)
            return 1

        name = entry["name"]
        group = entry.get("group")
        style_key = group or name
        style_label = group or entry["label"]
        pdg = entry["pdg"]
        max_as = entry["max_as"]
    
        if style_key not in grid_colors:
            style_index = len(grid_colors)
            grid_colors[style_key] = colors[style_index % len(colors)]
            grid_linestyles[style_key] = line_styles[style_index % len(line_styles)]
            legend_handles[style_key] = Patch(
                facecolor=grid_colors[style_key],
                edgecolor=grid_colors[style_key],
                linewidth=1.0,
                linestyle=grid_linestyles[style_key],
                label=style_label,
                alpha=0.55,
            )

        pdg_groups.setdefault(pdg, []).append(entry)

    ordered_pdgs = sorted(pdg_groups, key=lambda pdg: (-1 if pdg is None else pdg))
    ncols = len(ordered_pdgs)
    plot_width_in = 0.46 * 426.79135 / 72.27 
    fig, axes = plt.subplots(
        1,
        ncols,
        figsize=(plot_width_in, plot_width_in),
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes[0]
    heatmap_colorbars = {ax: [] for ax in axes}

    for ax, pdg in zip(axes, ordered_pdgs):
        grouped_support = {}

        for entry in pdg_groups[pdg]:
            aux_grid = entry["aux_grid"]
            lx_edges = np.asarray(aux_grid["lx_edges"], dtype=float)
            lq_edges = np.asarray(aux_grid["lq_edges"], dtype=float)
            weights = np.asarray(aux_grid["weights"], dtype=float)
            style_key = entry.get("group") or entry["name"]

            x_edges = np.power(10.0, lx_edges)
            q2_edges = np.power(10.0, lq_edges)

            if style_key not in grouped_support:
                grouped_support[style_key] = {
                    "x_edges": x_edges,
                    "q2_edges": q2_edges,
                    "selected": weights > 0.0,
                    "weights": weights.copy(),
                }
            else:
                grouped_support[style_key]["selected"] |= weights > 0.0
                grouped_support[style_key]["weights"] += weights

        # Draw the complete auxiliary grid independently of the selected cells.
        # It remains visible as orientation in the background, while the white
        # mesh edges below keep overlapping process regions easy to distinguish.
        reference_support = next(iter(grouped_support.values()))
        grid_x_edges = reference_support["x_edges"]
        grid_q2_edges = reference_support["q2_edges"]
        ax.vlines(
            grid_x_edges,
            grid_q2_edges[0],
            grid_q2_edges[-1],
            colors=(0.15, 0.15, 0.15, 0.28),
            linewidth=0.45,
            zorder=0,
        )
        ax.hlines(
            grid_q2_edges,
            grid_x_edges[0],
            grid_x_edges[-1],
            colors=(0.15, 0.15, 0.15, 0.28),
            linewidth=0.45,
            zorder=0,
        )

        for style_key, support in grouped_support.items():
            x_edges = support["x_edges"]
            q2_edges = support["q2_edges"]
            selected = support["selected"]
            weights = support["weights"]
            color = grid_colors[style_key]

            if plot_style == "boundary":
                mask = np.where(selected, 1.0, np.nan)
                ax.pcolormesh(
                    x_edges,
                    q2_edges,
                    mask,
                    shading="flat",
                    cmap=plt.matplotlib.colors.ListedColormap([color]),
                    vmin=0.0,
                    vmax=1.0,
                    # White cell edges keep the auxiliary-grid resolution visible.
                    edgecolors=(1.0, 1.0, 1.0, 0.18),
                    linewidth=0.25,
                    alpha=0.55,
                )
            else:
                visible_weights = np.ma.masked_where(~selected, weights)
                maximum = float(visible_weights.max())
                gradient = single_color_map(color, style_key)
                mesh = ax.pcolormesh(
                    x_edges,
                    q2_edges,
                    visible_weights,
                    shading="flat",
                    cmap=gradient,
                    vmin=0.0,
                    vmax=maximum,
                    #vmax=0.14,
                    edgecolors=(1.0, 1.0, 1.0, 0.18),
                    linewidth=0.25,
                    alpha=0.95,
                )
                heatmap_colorbars[ax].append(mesh)

            boundary = support_boundary_segments(
                x_edges,
                q2_edges,
                selected,
            )
            if boundary:
                ax.add_collection(
                    LineCollection(
                        boundary,
                        colors=[color],
                        linewidths=1.2,
                        linestyles=grid_linestyles[style_key],
                        alpha=0.95,
                        zorder=3,
                    )
                )

        draw_overlap(ax, grouped_support)

        ax.set_xscale("log")
        ax.set_yscale("log")
        title = pdg_title(pdg)
        if title is not None:
            ax.set_title(title)
        ax.set_xlabel(r"$x$")
        ax.tick_params(axis="both", which="major", width=1.0, length=6)
        ax.minorticks_off()

        for spine in ax.spines.values():
            spine.set_linewidth(1.1)

    axes[0].set_ylabel(r"$\mu^2\;[\mathrm{GeV}^2]$")

    title_as = as_title(max_as)
    figure_title = f"Kinematic Support — {title_as}"
    if subtitle:
        figure_title += "\n"
    #fig.suptitle(figure_title)
    if subtitle:
        fig.text(0.5, 0.94, subtitle, ha="center", va="top")

    axes[-1].legend(
        handles=list(legend_handles.values()),
        loc="best",
        frameon=False,
        labelspacing=0.35,
        handlelength=1.5,
        handletextpad=0.4,
        borderaxespad=0.3,
    )

    if plot_style == "heatmap":
        # Preserve the boundary-mode axes size; colorbars add physical width.
        fig.canvas.draw()
        positions = {ax: ax.get_position().frozen() for ax in axes}
        original_width, height = fig.get_size_inches()
        colorbar_slot_in = 0.65
        expanded_width = original_width + colorbar_slot_in * sum(
            len(meshes) for meshes in heatmap_colorbars.values()
        )
        fig.set_layout_engine(None)
        fig.set_size_inches(expanded_width, height)
        offset_in = 0.0
        for ax in axes:
            position = positions[ax]
            left_in = position.x0 * original_width + offset_in
            width_in = position.width * original_width
            ax.set_position([
                left_in / expanded_width, position.y0,
                width_in / expanded_width, position.height,
            ])
            for index, mesh in enumerate(heatmap_colorbars[ax]):
                colorbar_left_in = (
                    left_in + width_in + 0.08 + index * colorbar_slot_in
                )
                colorbar_ax = fig.add_axes([
                    colorbar_left_in / expanded_width, position.y0,
                    0.08 / expanded_width, position.height,
                ])
                colorbar = fig.colorbar(mesh, cax=colorbar_ax)
                colorbar.set_label(r"$W_z$" if mode == "supp" else r"$S_z$")
            offset_in += len(heatmap_colorbars[ax]) * colorbar_slot_in

    fig.savefig(output_pdf, bbox_inches="tight", pad_inches=0.022)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
