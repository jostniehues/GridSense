//! Helpers to transform extracted PineAPPL contributions into plot-ready data.

use std::collections::BTreeMap;
use serde::Serialize;
use thiserror::Error;
use pineappl::grid::KinematicsContribution;

const LEG_1: u8 = 0b01;
const LEG_2: u8 = 0b10;
const BOTH_LEGS: u8 = LEG_1 | LEG_2;
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
    #[error("zero support in selected bin {0} (A_B = 0)")]
    EmptyBinSupport(usize),
    #[error("support_leg must be 1 or 2, got {0}")]
    InvalidSupportLeg(u8),
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

impl AuxGrid {
    /// Add the cell weights of another auxiliary grid with the same binning.
    pub fn merge(&mut self, other: Self) {
        for (target_row, source_row) in self.weights.iter_mut().zip(other.weights) {
            for (target_weight, source_weight) in target_row.iter_mut().zip(source_row) {
                *target_weight += source_weight;
            }
        }
    }
}

pub(crate) struct BinSupportMaps {
    pub(crate) lx_edges: Vec<f64>,
    pub(crate) lq_edges: Vec<f64>,
    pub(crate) weights: Vec<Vec<Vec<f64>>>,
}

/// Signed contributions C_{B,u} for one observable bin.
/// Index order: weights[iq][ix1][ix2]
struct AuxGrid3D {
    weights: Vec<Vec<Vec<[f64; 4]>>>,
}

impl AuxGrid3D {
    fn new(nx: usize, nq: usize) -> Self {
        Self {
            weights: vec![vec![vec![[0.0; 4]; nx]; nx]; nq],
        }
    }

    fn add(&mut self, iq: usize, ix1: usize, ix2: usize, leg_mask: u8, value: f64) {
        self.weights[iq][ix1][ix2][usize::from(leg_mask)] += value;
    }

    fn absolute_sum(&self, requested_mask: u8) -> f64 {
        self.weights
            .iter()
            .flatten()
            .flatten()
            .flat_map(|cell| {
                [LEG_1, LEG_2, BOTH_LEGS].into_iter()
                    .filter(move |mask| mask & requested_mask != 0)
                    .map(|mask| cell[usize::from(mask)].abs())
            })
            .sum()
    }
}

fn logarithmic_edges(min: f64, max: f64, number_of_bins: usize) -> Vec<f64> {
    let log_min = min.log10();
    let log_max = max.log10();

    (0..=number_of_bins)
        .map(|i| { log_min + ((log_max - log_min) * (i as f64)) / (number_of_bins as f64) })
        .collect()
}

/// Sort the kinematic contribution from a given bin in the auxilary 3D-grid cell
fn accumulate_bin_in_3d(
    bin: &[KinematicsContribution],
    lx_edges: &[f64],
    lq_edges: &[f64]
) -> AuxGrid3D {
    let nx = lx_edges.len() - 1;
    let nq = lq_edges.len() - 1;

    let mut cells = AuxGrid3D::new(nx, nq);

    for point in bin {
        let Some(ix1) = find_bin(point.x1.log10(), lx_edges) else {
            continue;
        };

        let Some(ix2) = find_bin(point.x2.log10(), lx_edges) else {
            continue;
        };

        let Some(iq) = find_bin(point.q2.log10(), lq_edges) else {
            continue;
        };

        cells.add(iq, ix1, ix2, point.leg_mask, point.value);
    }

    cells
}

/// Implementation of P_B
fn project_bin_to_2d(cells: &AuxGrid3D, requested_mask: u8) -> Option<Vec<Vec<f64>>> {
    // Normalize only mask groups that contribute to the selected projection.
    // Keep the original groups separate when taking absolute values.
    let a_b = cells.absolute_sum(requested_mask);

    if a_b == 0.0 {
        return None;
    }

    let nq = cells.weights.len();
    let nx = cells.weights.first()?.len();

    let mut projected = vec![vec![0.0; nx]; nq];

    for iq in 0..nq {
        for ix1 in 0..nx {
            for ix2 in 0..nx {
                for leg_mask in [LEG_1, LEG_2, BOTH_LEGS] {
                    let effective_mask = leg_mask & requested_mask;
                    if effective_mask == 0 {
                        continue;
                    }

                    let value = cells.weights[iq][ix1][ix2][usize::from(leg_mask)].abs();

                    let number_of_legs = effective_mask.count_ones() as f64;
                    let normalized = value / (a_b * number_of_legs);

                    if (effective_mask & LEG_1) != 0 {
                        projected[iq][ix1] += normalized;
                    }

                    if (effective_mask & LEG_2) != 0 {
                        projected[iq][ix2] += normalized;
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
    support_leg: Option<u8>,
) -> Result<BinSupportMaps, TransformError> {
    let requested_mask = match support_leg {
        None => BOTH_LEGS,
        Some(1) => LEG_1,
        Some(2) => LEG_2,
        Some(leg) => return Err(TransformError::InvalidSupportLeg(leg)),
    };

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
    let mut weights = Vec::with_capacity(contrib.len());

    for (bin_index, bin) in contrib.iter().enumerate() {
        let cells = accumulate_bin_in_3d(bin, &lx_edges, &lq_edges);

        let projected = project_bin_to_2d(&cells, requested_mask)
            .ok_or(TransformError::EmptyBinSupport(bin_index))?;

        weights.push(projected);
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
pub fn average_support_maps(bin_support: BinSupportMaps) -> AuxGrid {
    let nx = bin_support.lx_edges.len() - 1;
    let nq = bin_support.lq_edges.len() - 1;
    let mut combined = vec![
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

#[cfg(test)]
mod tests {
    use super::*;

    fn spec() -> GridSpec {
        GridSpec {
            nx_x: 2, nx_q2: 1,
            x_min: 0.01, x_max: 1.0,
            q2_min: 1.0, q2_max: 100.0,
        }
    }

    fn point(x1: f64, x2: f64, value: f64, leg_mask: u8) -> KinematicsContribution {
        KinematicsContribution { x1, x2, q2: 10.0, value, leg_mask }
    }

    fn assert_row(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
        }
    }

    #[test]
    fn both_legs_reproduce_the_mean_of_normalized_single_leg_maps() {
        let contributions = vec![vec![
            point(0.02, 0.2, 3.0, BOTH_LEGS),
            point(0.02, 0.2, -1.0, BOTH_LEGS),
            point(0.02, 0.02, 1.0, BOTH_LEGS),
        ]];
        let both = build_per_bin_support_maps(&contributions, &spec(), None).unwrap();
        let first = build_per_bin_support_maps(&contributions, &spec(), Some(1)).unwrap();
        let second = build_per_bin_support_maps(&contributions, &spec(), Some(2)).unwrap();
        assert_row(&both.weights[0][0], &[2.0 / 3.0, 1.0 / 3.0]);
        assert_row(&first.weights[0][0], &[1.0, 0.0]);
        assert_row(&second.weights[0][0], &[1.0 / 3.0, 2.0 / 3.0]);
        for ix in 0..2 {
            assert!((both.weights[0][0][ix]
                - (first.weights[0][0][ix] + second.weights[0][0][ix]) / 2.0).abs() < 1e-12);
        }
    }

    #[test]
    fn flavour_masks_are_intersected_and_not_merged_before_absolute_values() {
        let contributions = vec![vec![
            point(0.02, 0.2, 100.0, LEG_1), // Excluded from leg 2 and its denominator.
            point(0.02, 0.2, 3.0, LEG_2),
            point(0.02, 0.2, -2.0, BOTH_LEGS), // Must not cancel the LEG_2 group.
            point(0.2, 0.02, 1.0, LEG_2),
        ]];
        let second = build_per_bin_support_maps(&contributions, &spec(), Some(2)).unwrap();
        assert_row(&second.weights[0][0], &[1.0 / 6.0, 5.0 / 6.0]);
    }

    #[test]
    fn zero_support_aborts_instead_of_dropping_a_bin() {
        for (empty_bin, leg) in [
            (vec![], None),
            (vec![point(0.02, 0.2, 1.0, LEG_1)], Some(2)),
            (vec![point(0.02, 0.2, 1.0, BOTH_LEGS),
                  point(0.02, 0.2, -1.0, BOTH_LEGS)], None),
        ] {
            let contributions = vec![vec![point(0.02, 0.2, 1.0, BOTH_LEGS)], empty_bin];
            match build_per_bin_support_maps(&contributions, &spec(), leg) {
                Err(TransformError::EmptyBinSupport(1)) => {},
                _ => panic!("expected zero-support error for second bin"),
            }
        }
    }
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
