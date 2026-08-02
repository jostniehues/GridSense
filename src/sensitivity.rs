//! Experimental inputs and covariance matrices for sensitivity calculations.

use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;

use anyhow::{ Context, Result, bail };
use serde::Deserialize;
use nalgebra::DMatrix;

use crate::transform::{ AuxGrid, BinSupportMaps };

/// Experimental input for one observable.
#[derive(Debug, Deserialize)]
pub(crate) struct ExperimentalData {
    observable: Observable,
    bins: Vec<ExperimentalBin>,
    uncertainties: Vec<UncertaintySource>,
    #[serde(rename = "grid_to_dat_norm")]
    grid_to_data_norm: Option<GridToDataNorm>,
}

#[derive(Debug, Deserialize)]
struct Observable {
    name: String,
    unit: String,
}

/// One bin of a one- or multi-differential distribution.
///
/// A one-dimensional bin has one entry in `ranges`; a two-dimensional bin
/// has two. The covariance construction is independent of that number.
#[derive(Debug, Deserialize)]
struct ExperimentalBin {
    id: usize,
    ranges: BTreeMap<String, BinRange>,
    value: f64,
}

#[derive(Debug, Deserialize)]
struct BinRange {
    low: f64,
    high: f64,
}

#[derive(Debug, Deserialize)]
struct UncertaintySource {
    name: String,
    kind: UncertaintyKind,
    unit: UncertaintyUnit,
    per_bin: Vec<f64>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum UncertaintyKind {
    Uncorrelated,
    Correlated,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum UncertaintyUnit {
    Absolute,
    Relative,
}

#[derive(Debug, Deserialize)]
struct GridToDataNorm {
    kind: String,
    factor: f64,
}

type CovarianceMatrix = DMatrix<f64>;

/// Read the YAML format (simplified HepData version).
pub(crate) fn read_experimental_data(path: &Path) -> Result<ExperimentalData> {
    let file = File::open(path).with_context(||
        format!("could not open experimental data '{}'", path.display())
    )?;

    serde_yaml
        ::from_reader(file)
        .with_context(|| format!("could not parse experimental data '{}'", path.display()))
}

/// Construct the covariance matrix from the uncertainty sources.
///
/// The first approximation accepts uncorrelated uncertainties and constructs
/// the covariance of relative deviations. Absolute uncertainties are first
/// divided by the experimental central value:
/// `V_ii = sum_s (delta_(s,i) / D_i)^2`, `V_ij = 0` for `i != j`.
/// Sources declared as relative are assumed to contain fractional values.
fn construct_covariance(data: &ExperimentalData) -> Result<CovarianceMatrix> {
    let number_of_bins = data.bins.len();
    let mut covariance = DMatrix::zeros(number_of_bins, number_of_bins);

    for source in &data.uncertainties {
        if source.per_bin.len() != number_of_bins {
            bail!(
                "uncertainty source '{}' has {} entries for {} bins",
                source.name,
                source.per_bin.len(),
                number_of_bins
            );
        }

        if source.kind != UncertaintyKind::Uncorrelated {
            bail!("uncertainty source '{}' is correlated, which is not supported yet", source.name);
        }

        for (bin_index, uncertainty) in source.per_bin.iter().enumerate() {
            if !uncertainty.is_finite() || *uncertainty < 0.0 {
                bail!("uncertainty source '{}' contains an invalid value", source.name);
            }

            let relative_uncertainty = match source.unit {
                UncertaintyUnit::Absolute => {
                    let central_value = data.bins[bin_index].value;

                    if !central_value.is_finite() || central_value == 0.0 {
                        bail!(
                            "cannot convert uncertainty source '{}' to relative units because an experimental central value is zero or non-finite",
                            source.name
                        );
                    }

                    uncertainty / central_value.abs()
                }
                UncertaintyUnit::Relative => *uncertainty,
            };

            covariance[(bin_index, bin_index)] += relative_uncertainty.powi(2);
        }
    }

    Ok(covariance)
}




fn invert_covariance(
    covariance: &CovarianceMatrix,
) -> Result<CovarianceMatrix> {
    if covariance.nrows() != covariance.ncols() {
        anyhow::bail!("covariance matrix is not square");
    }

    covariance
        .clone()
        .try_inverse()
        .ok_or_else(|| anyhow::anyhow!(
            "covariance matrix is singular and cannot be inverted"
        ))
}

fn build_p_matrix(
    support: &BinSupportMaps,
    number_of_data_bins: usize,
) -> Result<DMatrix<f64>> {
    if support.weights.len() != number_of_data_bins {
        bail!(
            "grid has {} observable bins but experimental data has {} bins",
            support.weights.len(),
            number_of_data_bins
        );
    }

    let nq = support.lq_edges.len() - 1;
    let nx = support.lx_edges.len() - 1;
    let mut p_matrix = DMatrix::zeros(number_of_data_bins, nq * nx);

    for observable_bin in 0..number_of_data_bins {
        for iq in 0..nq {
            for ix in 0..nx {
                let cell_index = iq * nx + ix;
                p_matrix[(observable_bin, cell_index)] =
                    support.weights[observable_bin][iq][ix];
            }
        }
    }

    Ok(p_matrix)
}

fn build_sensitivity_map(
    support: BinSupportMaps,
    inverse_covariance: &DMatrix<f64>,
) -> Result<AuxGrid> {
    if inverse_covariance.nrows() != inverse_covariance.ncols() {
        bail!("inverse covariance matrix is not square");
    }

    let p_matrix = build_p_matrix(&support, inverse_covariance.nrows())?;
    let nx = support.lx_edges.len() - 1;
    let nq = support.lq_edges.len() - 1;
    let mut weights = vec![vec![0.0; nx]; nq];

    for cell_index in 0..p_matrix.ncols() {
        let p_z = p_matrix.column(cell_index);
        let weighted_p_z = inverse_covariance * p_z;
        let sensitivity = p_z.dot(&weighted_p_z);
        let iq = cell_index / nx;
        let ix = cell_index % nx;

        weights[iq][ix] = sensitivity;
    }

    Ok(AuxGrid {
        lx_edges: support.lx_edges,
        lq_edges: support.lq_edges,
        weights,
    })
}

/// Given 
pub(crate) fn calculate_sensitivity_map(
    support: BinSupportMaps,
    data: &ExperimentalData,
) -> Result<AuxGrid> {
    let covariance = construct_covariance(data)?;
    let inverse_covariance = invert_covariance(&covariance)?;
    let mut sensitivity = build_sensitivity_map(support, &inverse_covariance)?;

    let total: f64 = sensitivity.weights.iter().flatten().sum();

    if !total.is_finite() || total <= 0.0 {
        bail!("sensitivity map has zero or invalid total weight");
    }
    //normalise S_z to 1
    for weight in sensitivity.weights.iter_mut().flatten() {
        *weight /= total;
    }

    Ok(sensitivity)
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructs_ttbar_table_21_covariance() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("input/ttbar/Table_21.yaml");
        let data = read_experimental_data(&path).expect("Table 21 should be readable");
        let covariance = construct_covariance(&data).expect("covariance should be constructible");

        println!("Covariance matrix:");

        for row in covariance.row_iter() {
            let formatted = row
                .iter()
                .map(|value| format!("{value:.6e}"))
                .collect::<Vec<_>>()
                .join(" ");

            println!("{formatted}");
        }

        assert_eq!(data.bins.len(), 10);
        assert_eq!(data.bins[0].ranges.len(), 1);
        assert_eq!(covariance.nrows(), 10);
        assert_eq!(covariance.ncols(), 10);

        let expected_first_diagonal =
            (0.00184_f64 / 0.0736).powi(2) + (0.00279_f64 / 0.0736).powi(2);
        assert!((covariance[(0, 0)] - expected_first_diagonal).abs() < 1.0e-15);
        assert_eq!(covariance[(0, 1)], 0.0);
        assert_eq!(covariance[(9, 0)], 0.0);
    }

    #[test]
    fn format_supports_two_dimensional_bins() {
        let yaml =
            r#"
observable:
  name: "d2sigma/dxdy"
  unit: "pb"
bins:
  - id: 0
    ranges:
      x: { low: 0.0, high: 1.0 }
      y: { low: 2.0, high: 3.0 }
    value: 4.0
uncertainties:
  - name: "stat"
    kind: "uncorrelated"
    unit: "absolute"
    per_bin: [0.5]
grid_to_dat_norm: null
"#;

        let data: ExperimentalData = serde_yaml::from_str(yaml).unwrap();
        let covariance = construct_covariance(&data).unwrap();

        assert_eq!(data.bins[0].ranges.len(), 2);
        assert_eq!(covariance, DMatrix::from_element(1, 1, 0.015625));
    }

    #[test]
    fn rejects_uncertainty_length_mismatch() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("input/ttbar/Table_21.yaml");
        let mut data = read_experimental_data(&path).unwrap();
        data.uncertainties[0].per_bin.pop();

        let error = construct_covariance(&data).unwrap_err();

        assert!(error.to_string().contains("9 entries for 10 bins"));
    }

    #[test]
    fn inverts_ttbar_covariance() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("input/ttbar/Table_21.yaml");
        let data = read_experimental_data(&path).unwrap();
        let covariance = construct_covariance(&data).unwrap();

        let inverse = invert_covariance(&covariance).unwrap();
        let identity = &covariance * inverse;

        assert!(identity.is_identity(1.0e-12));
    }

    #[test]
    fn builds_p_vectors_as_matrix_columns() {
        let support = BinSupportMaps {
            lx_edges: vec![-2.0, -1.0, 0.0],
            lq_edges: vec![2.0, 3.0],
            weights: vec![
                vec![vec![0.25, 0.75]],
                vec![vec![0.60, 0.40]],
            ],
        };

        let p_matrix = build_p_matrix(&support, 2).unwrap();

        assert_eq!(p_matrix.nrows(), 2);
        assert_eq!(p_matrix.ncols(), 2);
        assert_eq!(p_matrix.column(0).into_owned(), DMatrix::from_column_slice(2, 1, &[0.25, 0.60]));
        assert_eq!(p_matrix.column(1).into_owned(), DMatrix::from_column_slice(2, 1, &[0.75, 0.40]));
    }

    #[test]
    fn builds_sensitivity_aux_grid() {
        let support = BinSupportMaps {
            lx_edges: vec![-2.0, -1.0, 0.0],
            lq_edges: vec![2.0, 3.0],
            weights: vec![
                vec![vec![0.25, 0.75]],
                vec![vec![0.60, 0.40]],
            ],
        };
        let inverse_covariance = DMatrix::from_row_slice(
            2,
            2,
            &[4.0, 0.0, 0.0, 9.0],
        );

        let sensitivity = build_sensitivity_map(support, &inverse_covariance).unwrap();

        assert!((sensitivity.weights[0][0] - 3.49).abs() < 1.0e-12);
        assert!((sensitivity.weights[0][1] - 3.69).abs() < 1.0e-12);
        assert_eq!(sensitivity.lx_edges, vec![-2.0, -1.0, 0.0]);
        assert_eq!(sensitivity.lq_edges, vec![2.0, 3.0]);
    }
}
