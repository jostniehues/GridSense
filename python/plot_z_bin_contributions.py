"""Plot individual LHCb Z-bin contributions to the Z-only sensitivity."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from archiv.plot_lhcb_bin_contributions import (
    bin_contributions,
    construct_covariance,
    format_range,
    load_experimental_data,
    load_support,
)
from plot_z_sensitivity_selected_support import selected_bin_indices


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare equal-error and covariance-weighted contributions of "
            "individual LHCb Z rapidity bins to the Z-only sensitivity."
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
        default=PROJECT_ROOT / "output/plots/thesis/lhcb/z_bin_contributions.pdf",
    )
    parser.add_argument("--y-min", type=float, default=2.625)
    parser.add_argument("--y-max", type=float, default=3.375)
    parser.add_argument(
        "--diagonal-covariance",
        action="store_true",
        help="discard off-diagonal experimental covariance entries",
    )
    return parser.parse_args()


def make_plot(
    labels: list[str],
    equal: np.ndarray,
    experimental: np.ndarray,
    selected_bins: list[int],
    diagonal_covariance: bool,
    output_path: Path,
) -> None:
    positions = np.arange(len(labels), dtype=float)
    bar_width = 0.39
    selected_color = "#244A73"
    light_copper = "#F3D9C1"
    dark_copper = "#A64219"

    fig, ax = plt.subplots(figsize=(10.0, 4.8), constrained_layout=True)
    ax.axvspan(
        min(selected_bins) - 0.5,
        max(selected_bins) + 0.5,
        color=selected_color,
        alpha=0.075,
        linewidth=0.0,
        zorder=0,
    )
    ax.bar(
        positions - bar_width / 2,
        100.0 * equal,
        width=bar_width,
        color=light_copper,
        edgecolor=dark_copper,
        linewidth=1.05,
        label=r"Equal relative uncertainties, $C\propto I$",
        zorder=3,
    )
    ax.bar(
        positions + bar_width / 2,
        100.0 * experimental,
        width=bar_width,
        color=dark_copper,
        edgecolor=dark_copper,
        linewidth=0.9,
        label=(
            r"Experimental uncertainties, diagonal $C$"
            if diagonal_covariance
            else r"Experimental covariance, $C=C_{\mathrm{exp}}$"
        ),
        zorder=3,
    )

    ax.set_ylabel(
        r"Contribution to $Z\rightarrow\mu^+\mu^-$ $S_z$ [%]",
        fontsize=14,
    )
    ax.set_xlabel(r"$y_Z$", fontsize=14)
    ax.set_xticks(positions, labels)
    ax.tick_params(
        axis="x",
        labelrotation=45,
        labelsize=9,
        width=1.0,
        length=5,
    )
    ax.tick_params(axis="y", labelsize=11, width=1.0, length=5)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="0.82", linewidth=0.75)
    ax.set_xlim(-0.6, len(labels) - 0.4)
    ax.legend(loc="upper left", frameon=False, fontsize=10.5)

    selected_set = set(selected_bins)
    for index, tick_label in enumerate(ax.get_xticklabels()):
        tick_label.set_horizontalalignment("right")
        tick_label.set_rotation_mode("anchor")
        if index in selected_set:
            tick_label.set_color(selected_color)
            tick_label.set_fontweight("bold")

    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    arguments = parse_arguments()
    data = load_experimental_data(arguments.data)
    support = load_support(arguments.support, len(data["bins"]))
    covariance = construct_covariance(data)
    if arguments.diagonal_covariance:
        covariance = np.diag(np.diag(covariance))
    selected_bins = selected_bin_indices(data, arguments.y_min, arguments.y_max)

    equal = bin_contributions(support, np.eye(len(data["bins"])))
    experimental = bin_contributions(support, covariance)
    labels = [format_range(entry) for entry in data["bins"]]
    make_plot(
        labels,
        equal,
        experimental,
        selected_bins,
        arguments.diagonal_covariance,
        arguments.output,
    )

    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
