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

    fn absolute_sum(&self) -> f64 {
        self.weights
            .iter()
            .flatten()
            .flatten()
            .flatten()
            .map(|value| value.abs())
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
fn project_bin_to_2d(cells: &AuxGrid3D) -> Option<Vec<Vec<f64>>> {
    // A_B = sum_u |C_{B,u}|
    let a_b = cells.absolute_sum();

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
                    let value = cells.weights[iq][ix1][ix2][usize::from(leg_mask)].abs();

                    let number_of_legs = leg_mask.count_ones() as f64;
                    let normalized = value / (a_b * number_of_legs);

                    if (leg_mask & LEG_1) != 0 {
                        projected[iq][ix1] += normalized;
                    }

                    if (leg_mask & LEG_2) != 0 {
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
    grid_spec: &GridSpec
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
    let mut weights = Vec::with_capacity(contrib.len());

    for bin in contrib {
        let cells = accumulate_bin_in_3d(bin, &lx_edges, &lq_edges);

        let Some(projected) = project_bin_to_2d(&cells) else {
            continue;
        };

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
