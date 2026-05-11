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

/// Contribution relevance after per-bin normalization and summing over all bins.
#[derive(Clone, Debug)]
pub struct RelevancePoint {
    pub x: f64,
    pub q2: f64,
    pub relevance: f64,
}

#[derive(Debug, Error)]
pub enum TransformError {
    #[error("no contributions were available for x projection {0:?}")]
    EmptyProjection(XVar),
    #[error("no contributions were available to build relevance")]
    EmptyRelevanceInput,
    #[error("invalid auxiliary grid specification")]
    InvalidAuxgridSpec,
    #[error("no relevance points fell into the auxiliary grid")]
    EmptyAuxGrid,
    #[error("coverage threshold must lie in [0, 1], got {0}")]
    InvalidCoverageThreshold(f64),
    #[error("structure of x1 and x2 AuxGrid differs")]
    DifferentAuxGrids,
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
    xvar: XVar,
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


/// Builds a global relevance score for each `(x, q2, xvar)` point.
///
/// Within each bin, absolute contributions are normalized to the bin total.
/// The normalized weights of identical points are then summed over all bins.
pub fn build_relevance(
    contrib: &[Vec<SummedContribution>],
) -> Result<Vec<RelevancePoint>, TransformError> {
    let mut grouped = BTreeMap::<(u64, u64), f64>::new();

    for bin in contrib {
        // Per bin: for each node l we take abs(contrib).
        let total_per_bin: f64 = bin.iter().map(|point| point.value.abs()).sum();

        if total_per_bin == 0.0 {
            continue;
        }

        for point in bin {
            let key = (point.x.to_bits(), point.q2.to_bits());
            let normalized = point.value.abs() / total_per_bin;

            *grouped.entry(key).or_insert(0.0) += normalized;
        }
    }

    if grouped.is_empty() {
        return Err(TransformError::EmptyRelevanceInput);
    }

    Ok(grouped
        .into_iter()
        .map(|((x_bits, q2_bits), relevance)| RelevancePoint {
            x: f64::from_bits(x_bits),
            q2: f64::from_bits(q2_bits),
            relevance,
        })
        .collect())
} 

/// Specification of the plottet grid
pub struct GridSpec{
    pub nx_q2: usize, //number of nodes
    pub nx_x: usize,
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

pub fn sort_into_auxgrid(
    rel_point: &[RelevancePoint],
    grid_spec: &GridSpec,
) -> Result<AuxGrid, TransformError> {
    if grid_spec.nx_x == 0
        || grid_spec.nx_q2 == 0
        || grid_spec.x_min <= 0.0
        || grid_spec.q2_min <= 0.0
        || grid_spec.x_min >= grid_spec.x_max
        || grid_spec.q2_min >= grid_spec.q2_max
    {
        return Err(TransformError::InvalidAuxgridSpec);
    }

    let lx_min = grid_spec.x_min.log10();
    let lx_max = grid_spec.x_max.log10();
    let lq_min = grid_spec.q2_min.log10();
    let lq_max = grid_spec.q2_max.log10();

    let lx_edges: Vec<f64> = (0..=grid_spec.nx_x)
        .map(|i| lx_min + (lx_max - lx_min) * (i as f64) / (grid_spec.nx_x as f64))
        .collect();
    let lq_edges: Vec<f64> = (0..=grid_spec.nx_q2)
        .map(|i| lq_min + (lq_max - lq_min) * (i as f64) / (grid_spec.nx_q2 as f64))
        .collect();

    //2D auxiliary grid: rows in q2, columns in x.
    let mut weights = vec![vec![0.0; grid_spec.nx_x]; grid_spec.nx_q2];
    let mut filled = false;

    for point in rel_point {
        let x = point.x;
        let q2 = point.q2;
        let w = point.relevance;

        if !(x > 0.0
            && q2 > 0.0
            && x >= grid_spec.x_min
            && x <= grid_spec.x_max
            && q2 >= grid_spec.q2_min
            && q2 <= grid_spec.q2_max)
        {
            continue;
        }

        let lx = x.log10();
        let lq = q2.log10();

        let Some(ix) = find_bin(lx, &lx_edges) else {
            continue;
        };
        let Some(iy) = find_bin(lq, &lq_edges) else {
            continue;
        };

        weights[iy][ix] += w;
        filled = true;
    }

    if !filled {
        return Err(TransformError::EmptyAuxGrid);
    }

    Ok(AuxGrid {
        lx_edges,
        lq_edges,
        weights,
    })
}


fn find_bin(value: f64, edges: &[f64]) -> Option<usize> {
    if edges.len() < 2 || value < edges[0] || value > *edges.last()? {
        return None;
    }

    for i in 0..(edges.len() - 1) {
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
/// 
/// This is the function that gets also called in the main program
/// 
pub fn ranked_cells(mut aux_grid: AuxGrid, threshold: f64) -> Result<AuxGrid, TransformError> {
    if !(0.0..=1.0).contains(&threshold) {
        return Err(TransformError::InvalidCoverageThreshold(threshold));
    }

    let total_weight: f64 = aux_grid
        .weights
        .iter()
        .flatten()
        .copied()
        .sum();

    if total_weight == 0.0 {
        return Err(TransformError::EmptyAuxGrid);
    }

    let target_weight = threshold * total_weight;

    let mut ranked_cells: Vec<(usize, usize, f64)> = aux_grid
        .weights
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

/// Build the union of the x1 and x2 maps on the same auxiliary grid.
pub fn union_aux_grid(x1: &AuxGrid, x2: &AuxGrid) -> Result<AuxGrid, TransformError> {
    let lq_match = do_vecs_match(&x1.lq_edges, &x2.lq_edges);
    let lx_match = do_vecs_match(&x1.lx_edges, &x2.lx_edges);
    if !lx_match || !lq_match {
        return Err(TransformError::DifferentAuxGrids);
    }

    if x1.weights.len() != x2.weights.len() {
        return Err(TransformError::DifferentAuxGrids);
    }

    let mut weights = Vec::with_capacity(x1.weights.len());

    for (row_x1, row_x2) in x1.weights.iter().zip(&x2.weights) {
        if row_x1.len() != row_x2.len() {
            return Err(TransformError::DifferentAuxGrids);
        }

        let averaged_row = row_x1
            .iter()
            .zip(row_x2.iter())
            .map(|(w1, w2)| (w1 + w2) / 2.0)
            .collect();

        weights.push(averaged_row);
    }

    Ok(AuxGrid {
        lx_edges: x1.lx_edges.clone(),
        lq_edges: x1.lq_edges.clone(),
        weights,
    })
}

fn do_vecs_match<T: PartialEq>(a: &[T], b: &[T]) -> bool {
    a == b
}
