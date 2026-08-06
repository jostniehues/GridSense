//! Helpers to transform extracted PineAPPL contributions into plot-ready data.

use std::collections::BTreeMap;
use serde::Serialize;
use thiserror::Error;
use pineappl::grid::KinematicsContribution;

/// Contribution after summing over one of the momentum fractions.
#[derive(Clone, Debug)]
pub struct SummedContribution {
    pub x: f64,
    pub q2: f64,
    pub value: f64,
    pub xvar: XVar,
}

#[derive(Debug, Error)]
pub enum TransformError {
    #[error("no contributions were available for x projection {0:?}")] EmptyProjection(XVar),
    #[error("no contributions were available to build relevance")]
    EmptyRelevanceInput,
    #[error("invalid auxiliary grid specification")]
    InvalidAuxgridSpec,
    #[error("no relevance points fell into the auxiliary grid")]
    EmptyAuxGrid,
    #[error("coverage threshold must lie in [0, 1], got {0}")] InvalidCoverageThreshold(f64),
}

#[derive(Clone, Copy, Debug)]
pub enum XVar {
    X1,
    X2,
}

/// Sums the `(x1, x2, q2)` contributions per bin over one x-variable.
///
/// `XVar::X1` keeps `x1` and sums over `x2`.
/// `XVar::X2` keeps `x2` and sums over `x1`.
pub fn sum_over_x(
    contrib: &[Vec<KinematicsContribution>],
    xvar: XVar
) -> Result<Vec<Vec<SummedContribution>>, TransformError> {
    let summed_bins: Vec<Vec<SummedContribution>> = contrib
        .into_iter()
        .map(|bin| {
            let mut grouped = BTreeMap::<(u64, u64), f64>::new();

            for point in bin {
                let x_value = match xvar {
                    XVar::X1 => point.x1,
                    XVar::X2 => point.x2,
                };
                let key = (x_value.to_bits(), point.q2.to_bits());

                *grouped.entry(key).or_insert(0.0) += point.value;
            }

            grouped
                .into_iter()
                .map(|((x_bits, q2_bits), value)| SummedContribution {
                    x: f64::from_bits(x_bits),
                    q2: f64::from_bits(q2_bits),
                    value,
                    xvar,
                })
                .collect()
        })
        .collect();

    if summed_bins.iter().all(Vec::is_empty) {
        return Err(TransformError::EmptyProjection(xvar));
    }

    Ok(summed_bins)
}

/// Specification of the plottet grid
pub struct GridSpec {
    pub nx_q2: usize, //number of nodes
    pub nx_x: usize, // number of bins for both x1 and x2
    pub x_min: f64,
    pub x_max: f64,
    pub q2_min: f64,
    pub q2_max: f64,
}
#[derive(Debug, Serialize)]
pub struct AuxGrid {
    pub lx_edges: Vec<f64>,
    pub lq_edges: Vec<f64>,
    pub weights: Vec<Vec<f64>>,
}

pub(crate) struct BinSupportMaps {
    pub(crate) lx_edges: Vec<f64>,
    pub(crate) lq_edges: Vec<f64>,
    pub(crate) weights: Vec<Vec<Vec<f64>>>,
}

/// Signed contributions C_{B,u} for one observable bin.
/// Index order: weights[iq][ix1][ix2]
struct AuxGrid3D {
    weights: Vec<Vec<Vec<f64>>>,
}

impl AuxGrid3D {
    fn new(nx1: usize, nx2: usize, nq: usize) -> Self {
        Self {
            weights: vec![vec![vec![0.0; nx2]; nx1]; nq],
        }
    }

    fn add(&mut self, iq: usize, ix1: usize, ix2: usize, value: f64) {
        self.weights[iq][ix1][ix2] += value;
    }

    fn absolute_sum(&self) -> f64 {
        self.weights
            .iter()
            .flatten()
            .flatten()
            .map(|value| value.abs())
            .sum()
    }
}

/// One translated three-dimensional partition and its conservative mappings
/// back to the common reference grid.
struct ShiftedGrid {
    lx1_edges: Vec<f64>,
    lx2_edges: Vec<f64>,
    lq_edges: Vec<f64>,
    x1_to_reference: Vec<Vec<(usize, f64)>>,
    x2_to_reference: Vec<Vec<(usize, f64)>>,
    q_to_reference: Vec<Vec<(usize, f64)>>,
}

fn logarithmic_edges(min: f64, max: f64, number_of_bins: usize) -> Vec<f64> {
    let log_min = min.log10();
    let log_max = max.log10();

    (0..=number_of_bins)
        .map(|i| { log_min + ((log_max - log_min) * (i as f64)) / (number_of_bins as f64) })
        .collect()
}

/// Translate only the internal edges by a fraction of the uniform bin width.
/// The outer phase-space limits remain fixed.
fn shifted_internal_edges(reference_edges: &[f64], shift_fraction: f64) -> Vec<f64> {
    let mut shifted = reference_edges.to_vec();

    if reference_edges.len() < 3 || shift_fraction == 0.0 {
        return shifted;
    }

    let bin_width = reference_edges[1] - reference_edges[0];
    for edge in &mut shifted[1..reference_edges.len() - 1] {
        *edge += shift_fraction * bin_width;
    }

    shifted
}

/// For each source cell, calculate the fractions overlapping the reference
/// cells. The fractions of every source cell sum to one because both grids
/// have the same fixed outer boundaries.
fn overlap_fractions(
    source_edges: &[f64],
    reference_edges: &[f64],
) -> Vec<Vec<(usize, f64)>> {
    source_edges
        .windows(2)
        .map(|source| {
            let source_width = source[1] - source[0];

            reference_edges
                .windows(2)
                .enumerate()
                .filter_map(|(index, reference)| {
                    let overlap = source[1].min(reference[1])
                        - source[0].max(reference[0]);

                    (overlap > 0.0).then_some((index, overlap / source_width))
                })
                .collect()
        })
        .collect()
}

fn build_shifted_grids(
    lx_edges: &[f64],
    lq_edges: &[f64],
    translation_ensemble: bool,
) -> Vec<ShiftedGrid> {
    let shifts: &[f64] = if translation_ensemble {
        //&[-1.0 / 10.0, 0.0, 1.0 / 10.0]
        &[-1.0 / 3.0, 0.0, 1.0 / 3.0]
    } else {
        &[0.0]
    };
    let mut grids = Vec::with_capacity(shifts.len().pow(3));

    for &x1_shift in shifts {
        for &x2_shift in shifts {
            for &q_shift in shifts {
                let lx1_edges = shifted_internal_edges(lx_edges, x1_shift);
                let lx2_edges = shifted_internal_edges(lx_edges, x2_shift);
                let shifted_lq_edges = shifted_internal_edges(lq_edges, q_shift);

                grids.push(ShiftedGrid {
                    x1_to_reference: overlap_fractions(&lx1_edges, lx_edges),
                    x2_to_reference: overlap_fractions(&lx2_edges, lx_edges),
                    q_to_reference: overlap_fractions(&shifted_lq_edges, lq_edges),
                    lx1_edges,
                    lx2_edges,
                    lq_edges: shifted_lq_edges,
                });
            }
        }
    }

    grids
}

/// Sort the kinematic contribution from a given bin in the auxilary 3D-grid cell
fn accumulate_bin_in_3d(
    bin: &[KinematicsContribution],
    lx1_edges: &[f64],
    lx2_edges: &[f64],
    lq_edges: &[f64]
) -> AuxGrid3D {
    let nx1 = lx1_edges.len() - 1;
    let nx2 = lx2_edges.len() - 1;
    let nq = lq_edges.len() - 1;

    let mut cells = AuxGrid3D::new(nx1, nx2, nq);

    for point in bin {
        let Some(ix1) = find_bin(point.x1.log10(), lx1_edges) else {
            continue;
        };

        let Some(ix2) = find_bin(point.x2.log10(), lx2_edges) else {
            continue;
        };

        let Some(iq) = find_bin(point.q2.log10(), lq_edges) else {
            continue;
        };

        cells.add(iq, ix1, ix2, point.value);
    }

    cells
}

/// Build P_B on the common reference grid from one translated partition.
fn project_bin_to_reference_2d(
    cells: &AuxGrid3D,
    shifted_grid: &ShiftedGrid,
    reference_nx: usize,
    reference_nq: usize,
) -> Option<Vec<Vec<f64>>> {
    // A_B = sum_u |C_{B,u}|
    let a_b = cells.absolute_sum();

    if a_b == 0.0 {
        return None;
    }

    let nq = cells.weights.len();
    let nx1 = cells.weights.first()?.len();
    let nx2 = cells.weights.first()?.first()?.len();
    let mut projected = vec![vec![0.0; reference_nx]; reference_nq];

    for iq in 0..nq {
        for ix1 in 0..nx1 {
            for ix2 in 0..nx2 {
                let normalized = cells.weights[iq][ix1][ix2].abs() / (2.0 * a_b);

                for &(reference_iq, q_fraction) in &shifted_grid.q_to_reference[iq] {
                    for &(reference_ix, x_fraction) in
                        &shifted_grid.x1_to_reference[ix1]
                    {
                        projected[reference_iq][reference_ix] +=
                            normalized * q_fraction * x_fraction;
                    }

                    for &(reference_ix, x_fraction) in
                        &shifted_grid.x2_to_reference[ix2]
                    {
                        projected[reference_iq][reference_ix] +=
                            normalized * q_fraction * x_fraction;
                    }
                }
            }
        }
    }

    Some(projected)
}

/// Calculate all single P_B,z. 
pub(crate) fn build_per_bin_support_maps(
    contrib: &[Vec<KinematicsContribution>],
    grid_spec: &GridSpec,
    translation_ensemble: bool,
) -> Result<BinSupportMaps, TransformError> {
    if
        grid_spec.nx_x == 0 ||
        grid_spec.nx_q2 == 0 ||
        grid_spec.x_min <= 0.0 ||
        grid_spec.q2_min <= 0.0 ||
        grid_spec.x_min >= grid_spec.x_max ||
        grid_spec.q2_min >= grid_spec.q2_max
    {
        return Err(TransformError::InvalidAuxgridSpec);
    }

    let lx_edges = logarithmic_edges(grid_spec.x_min, grid_spec.x_max, grid_spec.nx_x);
    let lq_edges = logarithmic_edges(grid_spec.q2_min, grid_spec.q2_max, grid_spec.nx_q2);
    let shifted_grids = build_shifted_grids(
        &lx_edges,
        &lq_edges,
        translation_ensemble,
    );
    let mut weights = Vec::with_capacity(contrib.len());

    for bin in contrib {
        let mut averaged = vec![vec![0.0; grid_spec.nx_x]; grid_spec.nx_q2];
        let mut included_shifts = 0usize;

        for shifted_grid in &shifted_grids {
            let cells = accumulate_bin_in_3d(
                bin,
                &shifted_grid.lx1_edges,
                &shifted_grid.lx2_edges,
                &shifted_grid.lq_edges,
            );

            let Some(projected) = project_bin_to_reference_2d(
                &cells,
                shifted_grid,
                grid_spec.nx_x,
                grid_spec.nx_q2,
            ) else {
                continue;
            };

            included_shifts += 1;
            for iq in 0..grid_spec.nx_q2 {
                for ix in 0..grid_spec.nx_x {
                    averaged[iq][ix] += projected[iq][ix];
                }
            }
        }

        if included_shifts == 0 {
            continue;
        }

        for weight in averaged.iter_mut().flatten() {
            *weight /= included_shifts as f64;
        }

        weights.push(averaged);
    }

    if weights.is_empty() {
        return Err(TransformError::EmptyRelevanceInput);
    }

    Ok(BinSupportMaps {
        lx_edges,
        lq_edges,
        weights,
    })
}

/// Build W_Z = 1/N_B sum_B P_B,z
pub fn average_support_maps(
    bin_support: BinSupportMaps
) -> AuxGrid {
    let nx = bin_support.lx_edges.len() - 1;
    let nq = bin_support.lq_edges.len() - 1;
    let mut combined =
        vec![
            vec![0.0; nx];
            nq
        ];

    for projected in &bin_support.weights {
        for iq in 0..nq {
            for ix in 0..nx {
                combined[iq][ix] += projected[iq][ix];
            }
        }
    }

    for row in &mut combined {
        for weight in row {
            *weight /= bin_support.weights.len() as f64;
        }
    }

    AuxGrid {
        lx_edges: bin_support.lx_edges,
        lq_edges: bin_support.lq_edges,
        weights: combined,
    }
}

fn find_bin(value: f64, edges: &[f64]) -> Option<usize> {
    if edges.len() < 2 || value < edges[0] || value > *edges.last()? {
        return None;
    }

    for i in 0..edges.len() - 1 {
        if value >= edges[i] && value < edges[i + 1] {
            return Some(i);
        }
    }

    if value == *edges.last()? {
        return Some(edges.len() - 2);
    }

    None
}

/// Rank the grid cells per ser size and only take 90% of the total
pub fn ranked_cells(mut aux_grid: AuxGrid, threshold: f64) -> Result<AuxGrid, TransformError> {
    if !(0.0..=1.0).contains(&threshold) {
        return Err(TransformError::InvalidCoverageThreshold(threshold));
    }

    let total_weight: f64 = aux_grid.weights.iter().flatten().copied().sum();

    if total_weight == 0.0 {
        return Err(TransformError::EmptyAuxGrid);
    }

    let target_weight = threshold * total_weight;

    let mut ranked_cells: Vec<(usize, usize, f64)> = aux_grid.weights
        .iter()
        .enumerate()
        .flat_map(|(iy, row)| {
            row.iter()
                .enumerate()
                .map(move |(ix, &weight)| (iy, ix, weight))
        })
        .collect();

    // sort cell
    ranked_cells.sort_by(|a, b| b.2.total_cmp(&a.2));

    let mut cumulative = 0.0;

    // cumulate cells until threshold is reached
    for (rank, (iy, ix, weight)) in ranked_cells.into_iter().enumerate() {
        if rank > 0 && cumulative >= target_weight {
            aux_grid.weights[iy][ix] = 0.0;
            continue;
        }

        cumulative += weight;
    }

    Ok(aux_grid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translated_edges_keep_outer_boundaries_fixed() {
        let reference = vec![0.0, 1.0, 2.0, 3.0];
        let shifted = shifted_internal_edges(&reference, 1.0 / 3.0);

        assert_eq!(shifted[0], reference[0]);
        assert_eq!(shifted[3], reference[3]);
        assert!((shifted[1] - 4.0 / 3.0).abs() < 1.0e-12);
        assert!((shifted[2] - 7.0 / 3.0).abs() < 1.0e-12);
        assert!(shifted.windows(2).all(|edges| edges[0] < edges[1]));
    }

    #[test]
    fn overlap_projection_is_conservative() {
        let reference = vec![0.0, 1.0, 2.0, 3.0];
        let shifted = shifted_internal_edges(&reference, 1.0 / 3.0);
        let overlaps = overlap_fractions(&shifted, &reference);

        for source_cell in overlaps {
            let total: f64 = source_cell.iter().map(|(_, fraction)| fraction).sum();
            assert!((total - 1.0).abs() < 1.0e-12);
        }
    }

    #[test]
    fn ensemble_contains_all_three_dimensional_translations() {
        let x_edges = vec![-2.0, -1.0, 0.0];
        let q_edges = vec![2.0, 3.0, 4.0];

        assert_eq!(build_shifted_grids(&x_edges, &q_edges, false).len(), 1);
        assert_eq!(build_shifted_grids(&x_edges, &q_edges, true).len(), 27);
    }

    #[test]
    fn every_reprojected_support_map_remains_normalized() {
        let x_edges = vec![-2.0, -1.0, 0.0];
        let q_edges = vec![2.0, 3.0, 4.0];

        for shifted_grid in build_shifted_grids(&x_edges, &q_edges, true) {
            let mut cells = AuxGrid3D::new(2, 2, 2);
            cells.add(0, 0, 1, 2.0);
            cells.add(1, 1, 0, -1.0);

            let projected = project_bin_to_reference_2d(
                &cells,
                &shifted_grid,
                2,
                2,
            )
            .unwrap();
            let total: f64 = projected.iter().flatten().sum();

            assert!((total - 1.0).abs() < 1.0e-12);
            assert!(projected.iter().flatten().all(|weight| *weight >= 0.0));
        }
    }
}
