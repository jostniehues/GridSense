"""Compare full Z support and sensitivity with selected-bin support."""

from __future__ import annotations

import argparse
from pathlib import Path

from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from plot_auxgrid import support_boundary_segments
from archiv.plot_lhcb_bin_contributions import (
    construct_covariance,
    load_experimental_data,
    load_support_grid,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot the full covariance-weighted LHCb Z sensitivity and overlay "
            "the kinematic support of a selected y_Z interval."
        )
    )
    parser.add_argument(
        "--support",
        type=Path,
        default=PROJECT_ROOT / "output/data/support/z0_lhcb.csv",
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=PROJECT_ROOT / "input/wz/Z_lhcb.yaml",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            PROJECT_ROOT
            / "output/plots/thesis/lhcb/z_support_vs_sensitivity.pdf"
        ),
    )
    parser.add_argument("--y-min", type=float, default=2.625)
    parser.add_argument("--y-max", type=float, default=3.375)
    parser.add_argument("--threshold", type=float, default=0.90)
    parser.add_argument(
        "--diagonal-covariance",
        action="store_true",
        help="discard off-diagonal experimental covariance entries",
    )
    return parser.parse_args()


def ranked_mask(weights: np.ndarray, threshold: float) -> np.ndarray:
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must lie between zero and one")

    flat = weights.ravel()
    total = np.sum(flat)
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("map has a zero or invalid normalization")

    ranked = np.argsort(-flat, kind="stable")
    selected = np.zeros(flat.size, dtype=bool)
    cumulative = 0.0
    target = threshold * total

    for index in ranked:
        if np.any(selected) and cumulative >= target:
            break
        selected[index] = True
        cumulative += flat[index]

    return selected.reshape(weights.shape)


def selected_bin_indices(data: dict, y_min: float, y_max: float) -> list[int]:
    selected = []

    for index, entry in enumerate(data["bins"]):
        rapidity = entry.get("ranges", {}).get("y_Z")
        if rapidity is None:
            raise ValueError("every Z bin must contain a y_Z range")

        if rapidity["low"] >= y_min and rapidity["high"] <= y_max:
            selected.append(index)

    if not selected:
        raise ValueError(f"no complete y_Z bins lie in [{y_min}, {y_max}]")

    return selected


def calculate_maps(
    support: np.ndarray,
    covariance: np.ndarray,
    selected_bins: list[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    number_of_bins, number_of_q2_cells, number_of_x_cells = support.shape
    flattened = support.reshape(number_of_bins, -1)
    weighted = np.linalg.solve(covariance, flattened)
    sensitivity = np.sum(flattened * weighted, axis=0)
    sensitivity = sensitivity.reshape(number_of_q2_cells, number_of_x_cells)
    sensitivity /= np.sum(sensitivity)

    full_support = np.mean(support, axis=0)
    full_support /= np.sum(full_support)
    selected_support = np.mean(support[selected_bins], axis=0)
    selected_support /= np.sum(selected_support)
    return full_support, sensitivity, selected_support


def make_plot(
    full_support: np.ndarray,
    sensitivity: np.ndarray,
    selected_support: np.ndarray,
    log10_x_edges: np.ndarray,
    log10_q2_edges: np.ndarray,
    y_min: float,
    y_max: float,
    threshold: float,
    output_path: Path,
) -> None:
    full_support_mask = ranked_mask(full_support, threshold)
    sensitivity_mask = ranked_mask(sensitivity, threshold)
    support_mask = ranked_mask(selected_support, threshold)
    x_edges = np.power(10.0, log10_x_edges)
    q2_edges = np.power(10.0, log10_q2_edges)
    selected_color = "#244A73"
    heatmap_cmap = LinearSegmentedColormap.from_list(
        "contrasting_copper",
        ["#FFF8F1", "#DF8A46", "#8F2F0D"],
    )
    visible_full_support = np.ma.masked_where(~full_support_mask, full_support)
    visible_sensitivity = np.ma.masked_where(~sensitivity_mask, sensitivity)
    common_maximum = float(max(np.max(full_support), np.max(sensitivity)))

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.8, 5.8),
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )

    meshes = []
    for ax, weights in zip(
        axes,
        [visible_full_support, visible_sensitivity],
    ):
        ax.vlines(
            x_edges,
            q2_edges[0],
            q2_edges[-1],
            colors=(0.15, 0.15, 0.15, 0.22),
            linewidth=0.45,
            zorder=0,
        )
        ax.hlines(
            q2_edges,
            x_edges[0],
            x_edges[-1],
            colors=(0.15, 0.15, 0.15, 0.22),
            linewidth=0.45,
            zorder=0,
        )
        meshes.append(
            ax.pcolormesh(
                x_edges,
                q2_edges,
                weights,
                shading="flat",
                cmap=heatmap_cmap,
                vmin=0.0,
                vmax=common_maximum,
                edgecolors=(1.0, 1.0, 1.0, 0.18),
                linewidth=0.25,
                alpha=0.84,
                zorder=1,
            )
        )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$x$", fontsize=18)
        ax.tick_params(axis="both", which="major", labelsize=12.5, width=1.0, length=6)
        ax.minorticks_off()

        for spine in ax.spines.values():
            spine.set_linewidth(1.1)

    axes[0].set_ylabel(r"$\mu^2\;[\mathrm{GeV}^2]$", fontsize=18)
    axes[0].set_title(r"(a) Kinematic support, $W_z$", fontsize=16)
    axes[1].set_title(r"(b) Experimental sensitivity, $S_z$", fontsize=16)

    for ax, mesh, label in zip(axes, meshes, [r"$W_z$", r"$S_z$"]):
        colorbar = fig.colorbar(mesh, ax=ax, pad=0.02)
        colorbar.set_label(label, fontsize=14)
        colorbar.ax.tick_params(labelsize=10)

    boundary = support_boundary_segments(x_edges, q2_edges, support_mask)
    if boundary:
        axes[1].add_collection(
            LineCollection(
                boundary,
                colors=[selected_color],
                linewidths=2.5,
                linestyles="--",
                alpha=0.98,
                zorder=4,
            )
        )

    legend = [
        Line2D(
            [0],
            [0],
            color=selected_color,
            linewidth=2.5,
            linestyle="--",
            label=rf"Selected-bin $W_z$: ${y_min:g}<y_Z<{y_max:g}$",
        ),
    ]
    axes[1].legend(handles=legend, loc="best", frameon=False, fontsize=11.5)
    for ax in axes:
        ax.set_xlim(x_edges[0], x_edges[-1])
        ax.set_ylim(q2_edges[0], q2_edges[-1])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    arguments = parse_arguments()
    data = load_experimental_data(arguments.data)
    support_grid = load_support_grid(arguments.support, len(data["bins"]))
    covariance = construct_covariance(data)
    if arguments.diagonal_covariance:
        covariance = np.diag(np.diag(covariance))
    selected_bins = selected_bin_indices(data, arguments.y_min, arguments.y_max)
    full_support, sensitivity, selected_support = calculate_maps(
        support_grid.weights,
        covariance,
        selected_bins,
    )

    make_plot(
        full_support,
        sensitivity,
        selected_support,
        support_grid.log10_x_edges,
        support_grid.log10_q2_edges,
        arguments.y_min,
        arguments.y_max,
        arguments.threshold,
        arguments.output,
    )

    displayed_bins = ", ".join(str(index + 1) for index in selected_bins)
    print(f"wrote {arguments.output}")
    print(f"selected Z bins (one-based): {displayed_bins}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
