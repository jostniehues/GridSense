#!/usr/bin/env python3
"""Convert the CMS 13 TeV R=0.4 HEPData tables to GridSense YAML."""

from __future__ import annotations

import argparse
import tarfile
from pathlib import Path

import yaml


NUMBER_OF_RAPIDITY_BINS = 4


def read_yaml_member(archive: tarfile.TarFile, filename: str) -> dict:
    suffix = f"/{filename}"
    member = next(
        (entry for entry in archive.getmembers() if entry.name.endswith(suffix)),
        None,
    )
    if member is None:
        raise ValueError(f"archive does not contain {filename}")

    stream = archive.extractfile(member)
    if stream is None:
        raise ValueError(f"could not read {filename}")

    return yaml.safe_load(stream)


def format_number(value: float) -> str:
    return format(value, ".10g")


def load_measurement(archive: tarfile.TarFile):
    bins = []
    uncertainties = []
    systematic_shifts = {}
    bin_counts = []

    for rapidity_index in range(1, NUMBER_OF_RAPIDITY_BINS + 1):
        table = read_yaml_member(archive, f"ak4_xsec_ybin{rapidity_index}.yaml")
        independent = table["independent_variables"][0]["values"]
        dependent = table["dependent_variables"][0]["values"]

        if len(independent) != len(dependent):
            raise ValueError(f"xsec y-bin {rapidity_index} has inconsistent bin counts")

        bin_counts.append(len(dependent))
        uncertainty_label = f"all uncorr. unc. {rapidity_index}"

        for pt_bin, value in zip(independent, dependent):
            errors = {error["label"]: error["symerror"] for error in value["errors"]}
            uncertainty = next(
                (
                    shift
                    for label, shift in errors.items()
                    if label == uncertainty_label
                ),
                None,
            )
            if uncertainty is None:
                raise ValueError(
                    f"xsec y-bin {rapidity_index} is missing '{uncertainty_label}'"
                )

            bins.append((rapidity_index, pt_bin["low"], pt_bin["high"], value["value"]))
            uncertainties.append(abs(uncertainty))

            systematic_labels = {
                label
                for label in errors
                if not label.startswith("all uncorr. unc.") and label != "Lumi"
            }
            if systematic_shifts and systematic_labels != set(systematic_shifts):
                raise ValueError(
                    f"xsec y-bin {rapidity_index} has inconsistent systematic sources"
                )

            for label in systematic_labels:
                systematic_shifts.setdefault(label, []).append(errors[label])

    if any(len(shifts) != len(bins) for shifts in systematic_shifts.values()):
        raise ValueError("at least one systematic source does not cover every measured bin")

    return bins, uncertainties, systematic_shifts, bin_counts


def load_correlation_matrix(archive: tarfile.TarFile, bin_counts: list[int]):
    offsets = [0]
    for count in bin_counts:
        offsets.append(offsets[-1] + count)

    matrix = [[0.0] * offsets[-1] for _ in range(offsets[-1])]

    for row_group in range(NUMBER_OF_RAPIDITY_BINS):
        for column_group in range(NUMBER_OF_RAPIDITY_BINS):
            table = read_yaml_member(
                archive,
                f"ak4_correlations_y{row_group + 1}_y{column_group + 1}.yaml",
            )
            values = [
                item["value"] for item in table["dependent_variables"][0]["values"]
            ]
            rows = bin_counts[row_group]
            columns = bin_counts[column_group]

            if len(values) != rows * columns:
                raise ValueError(
                    f"correlation block y{row_group + 1}_y{column_group + 1} "
                    f"has {len(values)} entries, expected {rows * columns}"
                )

            for row in range(rows):
                for column in range(columns):
                    matrix[offsets[row_group] + row][offsets[column_group] + column] = (
                        values[row * columns + column]
                    )

    for row in range(len(matrix)):
        if abs(matrix[row][row] - 1.0) > 1.0e-12:
            raise ValueError(f"correlation matrix diagonal entry {row} is not one")

        for column in range(row):
            if abs(matrix[row][column] - matrix[column][row]) > 1.0e-12:
                raise ValueError(f"correlation matrix is not symmetric at ({row}, {column})")

    return matrix


def write_gridsense_yaml(
    output: Path,
    bins,
    uncertainties,
    systematic_shifts,
    bin_counts,
    correlation_matrix,
):
    lines = [
        "# CMS inclusive-jet cross section, anti-kT R=0.4, sqrt(s) = 13 TeV.",
        "# Source: HEPData record ins1972986, version 2, tables ak4_xsec_ybin1-4",
        "# and ak4_correlations_y1_y1 through ak4_correlations_y4_y4.",
        '# The published "all uncorr. unc." magnitudes are combined with the',
        "# particle-level correlation matrix from data and PYTHIA 8 unfolding.",
        "# All separately published systematic sources except Lumi are included",
        "# as fully correlated signed nuisance shifts.",
        "# Bin order: increasing |y| interval, then increasing pT within each interval.",
        "# The PineAPPL grids contain three additional unmeasured high-pT end bins.",
        f"grid_bin_counts: [{', '.join(map(str, bin_counts))}]",
        "",
        "bins:",
    ]

    for rapidity_index, low, high, value in bins:
        lines.append(
            f"  - value: {format_number(value)}  # |y| bin {rapidity_index}, "
            f"pT {format_number(low)}-{format_number(high)} GeV"
        )

    lines.extend(
        [
            "",
            "uncertainties:",
            '  - name: "data and unfolding"',
            '    kind: "correlated"',
            '    unit: "absolute"',
            "    per_bin: ["
            + ", ".join(format_number(value) for value in uncertainties)
            + "]",
            "    correlation_matrix:",
        ]
    )

    for row in correlation_matrix:
        lines.append("      - [" + ", ".join(format_number(value) for value in row) + "]")

    for name in sorted(systematic_shifts):
        lines.extend(
            [
                f'  - name: "{name}"',
                '    kind: "correlated"',
                '    unit: "absolute"',
                "    per_bin: ["
                + ", ".join(format_number(value) for value in systematic_shifts[name])
                + "]",
            ]
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument(
        "output",
        type=Path,
        help="output GridSense YAML, e.g. CMS_13TeV_R04_data_unfolding_correlated.yaml",
    )
    args = parser.parse_args()

    with tarfile.open(args.archive, "r:gz") as archive:
        bins, uncertainties, systematic_shifts, bin_counts = load_measurement(archive)
        correlation_matrix = load_correlation_matrix(archive, bin_counts)

    write_gridsense_yaml(
        args.output,
        bins,
        uncertainties,
        systematic_shifts,
        bin_counts,
        correlation_matrix,
    )


if __name__ == "__main__":
    main()
