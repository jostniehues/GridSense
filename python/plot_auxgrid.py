import json
import sys

from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np


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

    fig, ax = plt.subplots(figsize=(8.5, 6.8), constrained_layout=True)

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
    legend_handles = []

    for idx, entry in enumerate(data):
        if not isinstance(entry, dict):
            print("input JSON entry has unexpected structure", file=sys.stderr)
            return 1

        label = entry["label"]
        aux_grid = entry["aux_grid"]
        lx_edges = np.asarray(aux_grid["lx_edges"], dtype=float)
        lq_edges = np.asarray(aux_grid["lq_edges"], dtype=float)
        weights = np.asarray(aux_grid["weights"], dtype=float)

        x_edges = np.power(10.0, lx_edges)
        q2_edges = np.power(10.0, lq_edges)

        # Plot only the selected support area, not the raw relevance values.
        mask = np.where(weights > 0.0, 1.0, np.nan)

        ax.pcolormesh(
            x_edges,
            q2_edges,
            mask,
            shading="flat",
            cmap=plt.matplotlib.colors.ListedColormap([colors[idx % len(colors)]]),
            vmin=0.0,
            vmax=1.0,
            edgecolors=(0.0, 0.0, 0.0, 0.28),
            linewidth=0.22,
            alpha=0.65,
        )
        legend_handles.append(
            Patch(facecolor=colors[idx % len(colors)], edgecolor="none", label=label, alpha=0.9)
        )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_title("Grid-based Kinematic Support", fontsize=20)
    ax.set_xlabel(r"$x$", fontsize=18)
    ax.set_ylabel(r"$Q^2\;[\mathrm{GeV}^2]$", fontsize=18)

    # ax.grid(which="major", color="black", linewidth=0.8, alpha=0.45)
    # ax.grid(which="minor", color="black", linewidth=0.45, alpha=0.22)
    ax.tick_params(axis="both", which="major", labelsize=13, width=1.0, length=6)
    # ax.tick_params(axis="both", which="minor", width=0.8, length=3)

    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=True,
            fontsize=11,
            borderaxespad=0.0,
        )

    fig.savefig(output_pdf)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
