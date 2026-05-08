import json
import sys

from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np


def pdg_title(pdg):
    if pdg is None:
        return "All Channels"
    pdg_dic = {
    0: "Gluon" ,   
    1: "Down-Quark",
    2: "Up-Quark",
    3: "Strange-Quark",
    4: "Charm-Quark",
    5: "Bottom-Quark",
    6: "Top-Quark",

    -1: "Anti-Down-Quark",
    -2: "Anti-Up-Quark",
    -3: "Anti-Strange-Quark",
    -4: "Anti-Charm-Quark",
    -5: "Anti-Bottom-Quark",
    -6: "Anti-Top-Quark",
}

    return f"{pdg_dic[pdg]}"





def main() -> int:
    if len(sys.argv) != 3:
        print("usage: python/plot_auxgrid.py <input.json> <output.pdf>", file=sys.stderr)
        return 1

    input_json = sys.argv[1]
    output_pdf = sys.argv[2]

    with open(input_json, encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, list) or len(data) == 0:
        print("input JSON does not contain any auxiliary grids", file=sys.stderr)
        return 1

    colors = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#17becf",
    ]

    pdg_groups = {}
    grid_colors = {}
    legend_handles = {}

    for entry in data:
        if not isinstance(entry, dict):
            print("input JSON entry has unexpected structure", file=sys.stderr)
            return 1

        name = entry["name"]
        pdg = entry["pdg"]

        if name not in grid_colors:
            grid_colors[name] = colors[len(grid_colors) % len(colors)]
            legend_handles[name] = Patch(
                facecolor=grid_colors[name],
                edgecolor="none",
                label=entry["label"],
                alpha=0.9,
            )

        pdg_groups.setdefault(pdg, []).append(entry)

    ordered_pdgs = sorted(pdg_groups, key=lambda pdg: (-1 if pdg is None else pdg))
    ncols = len(ordered_pdgs)
    fig, axes = plt.subplots(
        1,
        ncols,
        figsize=(6.4 * ncols, 6.8),
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes[0]

    for ax, pdg in zip(axes, ordered_pdgs):
        for entry in pdg_groups[pdg]:
            aux_grid = entry["aux_grid"]
            lx_edges = np.asarray(aux_grid["lx_edges"], dtype=float)
            lq_edges = np.asarray(aux_grid["lq_edges"], dtype=float)
            weights = np.asarray(aux_grid["weights"], dtype=float)

            x_edges = np.power(10.0, lx_edges)
            q2_edges = np.power(10.0, lq_edges)

            # Plot only the selected support area, not the raw relevance values.
            mask = np.where(weights > 0.0, 1.0, np.nan)
            color = grid_colors[entry["name"]]

            ax.pcolormesh(
                x_edges,
                q2_edges,
                mask,
                shading="flat",
                cmap=plt.matplotlib.colors.ListedColormap([color]),
                vmin=0.0,
                vmax=1.0,
                edgecolors=(0.0, 0.0, 0.0, 0.28),
                linewidth=0.22,
                alpha=0.65,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(pdg_title(pdg), fontsize=16)
        ax.set_xlabel(r"$x$", fontsize=18)
        ax.tick_params(axis="both", which="major", labelsize=13, width=1.0, length=6)

        for spine in ax.spines.values():
            spine.set_linewidth(1.1)

    axes[0].set_ylabel(r"$Q^2\;[\mathrm{GeV}^2]$", fontsize=18)
    fig.suptitle("Grid-based Kinematic Support", fontsize=20)

    fig.legend(
        handles=list(legend_handles.values()),
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        bbox_transform=fig.transFigure,
        frameon=True,
        fontsize=14,
        borderaxespad=0.0,
    )

    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
