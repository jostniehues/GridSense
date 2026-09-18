"""Plot LHCb bin contributions to the grouped W/Z sensitivity map.

The script compares an equal-relative-uncertainty reference, C proportional
to the identity, with the experimental covariance assembled from the YAML
inputs.  It reads the unranked per-bin P_B,z maps exported by GridSense.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class InputSpec:
    name: str
    panel_title: str
    x_label: str
    data_path: Path
    support_path: Path


@dataclass(frozen=True)
class SupportGrid:
    weights: np.ndarray
    log10_x_edges: np.ndarray
    log10_q2_edges: np.ndarray


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare equal-error and covariance-weighted bin contributions "
            "to the grouped LHCb W/Z sensitivity."
        )
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "output/plots/thesis/lhcb/lhcb_bin_contributions.pdf",
        help="output figure (default: %(default)s)",
    )
    parser.add_argument(
        "--wm-support",
        type=Path,
        default=PROJECT_ROOT / "output/data/support/wm_lhcb.csv",
    )
    parser.add_argument(
        "--wp-support",
        type=Path,
        default=PROJECT_ROOT / "output/data/support/wp_lhcb.csv",
    )
    parser.add_argument(
        "--z-support",
        type=Path,
        default=PROJECT_ROOT / "output/data/support/z0_lhcb.csv",
    )
    parser.add_argument(
        "--wm-data",
        type=Path,
        default=PROJECT_ROOT / "input/wz/Wminus_lhcb.yaml",
    )
    parser.add_argument(
        "--wp-data",
        type=Path,
        default=PROJECT_ROOT / "input/wz/Wplus_lhcb.yaml",
    )
    parser.add_argument(
        "--z-data",
        type=Path,
        default=PROJECT_ROOT / "input/wz/Z_lhcb.yaml",
    )
    return parser.parse_args()


def load_experimental_data(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    if not isinstance(data, dict) or "bins" not in data or "uncertainties" not in data:
        raise ValueError(f"experimental input has an unexpected structure: {path}")

    return data


def load_support_grid(path: Path, number_of_bins: int) -> SupportGrid:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError(f"support file is empty: {path}")

    number_of_q2_cells = max(int(row["iq"]) for row in rows) + 1
    number_of_x_cells = max(int(row["ix"]) for row in rows) + 1
    support = np.zeros((number_of_bins, number_of_q2_cells, number_of_x_cells))
    log10_x_edges = np.full(number_of_x_cells + 1, np.nan)
    log10_q2_edges = np.full(number_of_q2_cells + 1, np.nan)

    for row in rows:
        observable_bin = int(row["observable_bin"])
        iq = int(row["iq"])
        ix = int(row["ix"])

        if observable_bin >= number_of_bins:
            raise ValueError(
                f"support file {path} contains observable bin {observable_bin}, "
                f"but the YAML contains only {number_of_bins} bins"
            )

        support[observable_bin, iq, ix] = float(row["p_bz"])
        log10_x_edges[ix] = float(row["log10_x_low"])
        log10_x_edges[ix + 1] = float(row["log10_x_high"])
        log10_q2_edges[iq] = float(row["log10_q2_low"])
        log10_q2_edges[iq + 1] = float(row["log10_q2_high"])

    if np.any(~np.isfinite(log10_x_edges)) or np.any(~np.isfinite(log10_q2_edges)):
        raise ValueError(f"support file has incomplete auxiliary-grid edges: {path}")

    return SupportGrid(support, log10_x_edges, log10_q2_edges)


def load_support(path: Path, number_of_bins: int) -> np.ndarray:
    support = load_support_grid(path, number_of_bins)
    return support.weights.reshape(number_of_bins, -1)


def construct_covariance(data: dict) -> np.ndarray:
    central_values = np.asarray([entry["value"] for entry in data["bins"]], dtype=float)
    number_of_bins = central_values.size
    covariance = np.zeros((number_of_bins, number_of_bins))

    for source in data["uncertainties"]:
        shifts = np.asarray(source["per_bin"], dtype=float)
        if shifts.size != number_of_bins:
            raise ValueError(
                f"uncertainty source {source['name']!r} has {shifts.size} values "
                f"for {number_of_bins} bins"
            )

        if source["unit"] == "absolute":
            shifts = shifts / np.abs(central_values)
        elif source["unit"] != "relative":
            raise ValueError(f"unknown uncertainty unit: {source['unit']!r}")

        correlation = source.get("correlation_matrix")
        if correlation is not None:
            correlation = np.asarray(correlation, dtype=float)
            if correlation.shape != (number_of_bins, number_of_bins):
                raise ValueError(
                    f"correlation matrix for {source['name']!r} has shape "
                    f"{correlation.shape}, expected {(number_of_bins, number_of_bins)}"
                )
            covariance += correlation * np.outer(shifts, shifts)
        elif source["kind"] == "uncorrelated":
            covariance += np.diag(shifts**2)
        elif source["kind"] == "correlated":
            covariance += np.outer(shifts, shifts)
        else:
            raise ValueError(f"unknown uncertainty kind: {source['kind']!r}")

    return covariance


def block_diagonal(matrices: list[np.ndarray]) -> np.ndarray:
    total_size = sum(matrix.shape[0] for matrix in matrices)
    result = np.zeros((total_size, total_size))
    offset = 0

    for matrix in matrices:
        size = matrix.shape[0]
        result[offset : offset + size, offset : offset + size] = matrix
        offset += size

    return result


def bin_contributions(
    support: np.ndarray,
    covariance: np.ndarray,
) -> np.ndarray:
    weighted_support = np.linalg.solve(covariance, support)
    attribution = np.sum(support * weighted_support, axis=1)
    total = np.sum(attribution)

    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("grouped sensitivity has a zero or invalid normalization")

    return attribution / total


def format_range(entry: dict) -> str:
    ranges = entry.get("ranges", {})
    if len(ranges) != 1:
        raise ValueError("the Z-bin plot expects exactly one range per observable bin")

    limits = next(iter(ranges.values()))
    return f"{limits['low']:g}–{limits['high']:g}"


def make_plot(
    panels: list[tuple[InputSpec, list[str], np.ndarray, np.ndarray]],
    output_path: Path,
) -> None:
    bar_width = 0.39
    color = "#665A8E"
    light_color = "#D9D5E5"

    fig, axes = plt.subplots(
        1,
        len(panels),
        figsize=(15.2, 5.2),
        constrained_layout=True,
        sharey=True,
    )

    legend_handles = None
    for ax, (specification, labels, equal, experimental) in zip(axes, panels):
        positions = np.arange(len(labels), dtype=float)
        ax.bar(
            positions - bar_width / 2,
            100.0 * equal,
            width=bar_width,
            color=light_color,
            edgecolor=color,
            linewidth=1.25,
            label=r"Equal relative uncertainties, $C\propto I$",
            zorder=3,
        )
        ax.bar(
            positions + bar_width / 2,
            100.0 * experimental,
            width=bar_width,
            color=color,
            edgecolor=color,
            linewidth=1.0,
            label=r"Experimental covariance, $C=C_{\mathrm{exp}}$",
            zorder=3,
        )

        ax.set_title(specification.panel_title, fontsize=16)
        ax.set_xlabel(specification.x_label, fontsize=14)
        ax.set_xticks(
            positions,
            labels,
            rotation=55,
            ha="right",
            rotation_mode="anchor",
        )
        ax.tick_params(axis="x", labelsize=7.5 if len(labels) > 10 else 9)
        ax.tick_params(axis="y", labelsize=11, width=1.0, length=5)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=(0.25, 0.25, 0.25, 0.14), linewidth=0.7)
        ax.axhline(0.0, color=(0.15, 0.15, 0.15, 0.7), linewidth=0.8)

        for spine in ax.spines.values():
            spine.set_linewidth(1.05)

        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()

    axes[0].set_ylabel(r"Contribution to grouped $S_z$ [\%]", fontsize=15)
    fig.suptitle(r"LHCb $W^\pm,Z$: bin contributions to grouped sensitivity", fontsize=17)
    fig.legend(
        legend_handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=2,
        frameon=False,
        fontsize=11.5,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    arguments = parse_arguments()
    specifications = [
        InputSpec("W-", r"$W^-$", r"$\eta_\mu$ bin", arguments.wm_data, arguments.wm_support),
        InputSpec("W+", r"$W^+$", r"$\eta_\mu$ bin", arguments.wp_data, arguments.wp_support),
        InputSpec("Z", r"$Z$", r"$y_Z$ bin", arguments.z_data, arguments.z_support),
    ]

    supports = []
    covariances = []
    data_sets = []

    for specification in specifications:
        data = load_experimental_data(specification.data_path)
        support = load_support(specification.support_path, len(data["bins"]))

        if supports and support.shape[1] != supports[0].shape[1]:
            raise ValueError("all support files must use the same auxiliary grid")

        data_sets.append(data)
        supports.append(support)
        covariances.append(construct_covariance(data))

    grouped_support = np.concatenate(supports, axis=0)
    experimental_covariance = block_diagonal(covariances)
    equal_covariance = np.eye(grouped_support.shape[0])

    equal = bin_contributions(grouped_support, equal_covariance)
    experimental = bin_contributions(grouped_support, experimental_covariance)

    contributions = {}
    offset = 0
    for specification, data in zip(specifications, data_sets):
        stop = offset + len(data["bins"])
        contributions[specification.name] = (
            specification,
            [format_range(entry) for entry in data["bins"]],
            equal[offset:stop],
            experimental[offset:stop],
        )
        offset = stop

    panels = [
        contributions["Z"],
        contributions["W+"],
        contributions["W-"],
    ]
    make_plot(panels, arguments.output)

    print(f"wrote {arguments.output}")
    for specification, labels, equal_process, experimental_process in panels:
        strongest = int(np.argmax(experimental_process - equal_process))
        print(
            f"{specification.name} share of grouped S_z: "
            f"equal = {100.0 * np.sum(equal_process):.2f}%, "
            f"experimental = {100.0 * np.sum(experimental_process):.2f}%; "
            f"largest increase in bin {labels[strongest]} "
            f"({100.0 * (experimental_process[strongest] - equal_process[strongest]):+.2f} pp)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
