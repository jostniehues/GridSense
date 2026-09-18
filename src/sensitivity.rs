//! Experimental inputs and covariance matrices for sensitivity calculations.

use std::collections::{ HashMap, HashSet };
use std::fs::File;
use std::path::Path;

use anyhow::{ Context, Result, bail };
use serde::Deserialize;
use nalgebra::DMatrix;

use crate::transform::{ AuxGrid, BinSupportMaps };

// Debug switch: treat every uncertainty contribution as bin-wise uncorrelated.
const IGNORE_CORRELATIONS: bool = false;

/// Experimental input for one observable.
#[derive(Debug, Deserialize)]
pub(crate) struct ExperimentalData {
    bins: Vec<ExperimentalBin>,
    uncertainties: Vec<UncertaintySource>,

    /// Optional number of measured bins retained from each configured grid
    /// Removes trailing theory-only bins before support maps are built
    #[serde(default)]
    grid_bin_counts: Option<Vec<usize>>,
}

impl ExperimentalData {
    pub(crate) fn grid_bin_counts(&self) -> Option<&[usize]> {
        self.grid_bin_counts.as_deref()
    }
}

/// One bin of a one- or multi-differential distribution.
///
/// A one-dimensional bin has one entry in `ranges`; a two-dimensional bin
/// has two. The covariance construction is independent of that number.
#[derive(Debug, Deserialize)]
struct ExperimentalBin {
    value: f64,
}

#[derive(Debug, Deserialize)]
struct UncertaintySource {
    name: String,
    kind: UncertaintyKind,
    unit: UncertaintyUnit,
    per_bin: Vec<f64>,
    #[serde(default)]
    correlation_matrix: Option<Vec<Vec<f64>>>,
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

/// Construct the covariance matrix of relative deviations.
/// Absolute shifts are divided by the experimental central values. Relative shifts are assumed to be fractional values.
fn construct_covariance(data: &ExperimentalData) -> Result<CovarianceMatrix> {
    construct_joint_covariance(&[data])
}

fn relative_shifts(data: &ExperimentalData, source: &UncertaintySource) -> Result<Vec<f64>> {
    let number_of_bins = data.bins.len();

    if source.per_bin.len() != number_of_bins {
        bail!(
            "uncertainty source '{}' has {} entries for {} bins",
            source.name,
            source.per_bin.len(),
            number_of_bins
        );
    }

    let mut shifts = Vec::with_capacity(number_of_bins);

    for (bin_index, uncertainty) in source.per_bin.iter().enumerate() {
        if !uncertainty.is_finite() {
            bail!("uncertainty source '{}' contains an invalid value", source.name);
        }

        if
            (source.kind == UncertaintyKind::Uncorrelated || source.correlation_matrix.is_some()) &&
            *uncertainty < 0.0
        {
            bail!("uncertainty source '{}' contains a negative uncertainty magnitude", source.name);
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

        shifts.push(relative_uncertainty);
    }

    Ok(shifts)
}

/// Construct one relative covariance matrix for several data sets. Correlated
/// sources with the same name share one nuisance parameter across data sets.
fn construct_joint_covariance(data_sets: &[&ExperimentalData]) -> Result<CovarianceMatrix> {
    let total_bins: usize = data_sets
        .iter()
        .map(|data| data.bins.len())
        .sum();
    let mut covariance = DMatrix::zeros(total_bins, total_bins);
    let mut shared_shifts = HashMap::<String, Vec<f64>>::new();
    let mut offset = 0;

    for data in data_sets {
        let number_of_bins = data.bins.len();
        let mut seen_shared_sources = HashSet::new();

        for source in &data.uncertainties {
            let shifts = relative_shifts(data, source)?;

            if IGNORE_CORRELATIONS {
                for (bin_index, uncertainty) in shifts.iter().enumerate() {
                    covariance[(offset + bin_index, offset + bin_index)] += uncertainty.powi(2);
                }
                continue;
            }

            if let Some(correlation) = &source.correlation_matrix {
                if source.kind != UncertaintyKind::Correlated {
                    bail!(
                        "uncertainty source '{}' has a correlation matrix but is not declared correlated",
                        source.name
                    );
                }

                validate_correlation_matrix(correlation, number_of_bins, &source.name)?;

                for row in 0..number_of_bins {
                    for column in 0..number_of_bins {
                        covariance[(offset + row, offset + column)] +=
                            correlation[row][column] * shifts[row] * shifts[column];
                    }
                }

                continue;
            }

            match source.kind {
                UncertaintyKind::Uncorrelated => {
                    for (bin_index, uncertainty) in shifts.iter().enumerate() {
                        covariance[(offset + bin_index, offset + bin_index)] += uncertainty.powi(2);
                    }
                }
                UncertaintyKind::Correlated => {
                    if !seen_shared_sources.insert(source.name.as_str()) {
                        bail!(
                            "correlated uncertainty source '{}' occurs more than once in one data set",
                            source.name
                        );
                    }

                    let combined = shared_shifts
                        .entry(source.name.clone())
                        .or_insert_with(|| vec![0.0; total_bins]);

                    for (bin_index, shift) in shifts.iter().enumerate() {
                        combined[offset + bin_index] = *shift;
                    }
                }
            }
        }

        offset += number_of_bins;
    }

    for shifts in shared_shifts.values() {
        for row in 0..total_bins {
            for column in 0..total_bins {
                covariance[(row, column)] += shifts[row] * shifts[column];
            }
        }
    }

    Ok(covariance)
}

fn validate_correlation_matrix(
    correlation: &[Vec<f64>],
    number_of_bins: usize,
    source_name: &str
) -> Result<()> {
    if
        correlation.len() != number_of_bins ||
        correlation.iter().any(|row| row.len() != number_of_bins)
    {
        bail!(
            "correlation matrix for source '{}' must be {}x{}",
            source_name,
            number_of_bins,
            number_of_bins
        );
    }

    for row in 0..number_of_bins {
        for column in 0..number_of_bins {
            let coefficient = correlation[row][column];

            if !coefficient.is_finite() || !(-1.0..=1.0).contains(&coefficient) {
                bail!("correlation matrix for source '{}' contains an invalid coefficient", source_name);
            }

            if (coefficient - correlation[column][row]).abs() > 1.0e-12 {
                bail!("correlation matrix for source '{}' is not symmetric", source_name);
            }
        }

        if (correlation[row][row] - 1.0).abs() > 1.0e-12 {
            bail!("correlation matrix for source '{}' must have ones on its diagonal", source_name);
        }
    }

    Ok(())
}

fn invert_covariance(covariance: &CovarianceMatrix) -> Result<CovarianceMatrix> {
    if covariance.nrows() != covariance.ncols() {
        anyhow::bail!("covariance matrix is not square");
    }

    covariance
        .clone()
        .try_inverse()
        .ok_or_else(|| anyhow::anyhow!("covariance matrix is singular and cannot be inverted"))
}

fn build_p_matrix(support: &BinSupportMaps, number_of_data_bins: usize) -> Result<DMatrix<f64>> {
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
                //flatten
                let cell_index = iq * nx + ix;
                p_matrix[(observable_bin, cell_index)] = support.weights[observable_bin][iq][ix];
            }
        }
    }

    Ok(p_matrix)
}

fn build_sensitivity_map(
    support: BinSupportMaps,
    inverse_covariance: &DMatrix<f64>
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
        //deflatten
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

pub(crate) fn calculate_sensitivity_map(
    support: BinSupportMaps,
    data: &ExperimentalData
) -> Result<AuxGrid> {
    let covariance = construct_covariance(data)?;
    let inverse_covariance = invert_covariance(&covariance)?;
    let sensitivity = build_sensitivity_map(support, &inverse_covariance)?;

    normalize_sensitivity(sensitivity)
}

pub(crate) fn calculate_joint_sensitivity_map(
    data_sets: Vec<(BinSupportMaps, ExperimentalData)>
) -> Result<AuxGrid> {
    if data_sets.is_empty() {
        bail!("cannot calculate sensitivity without data sets");
    }

    let experimental_data: Vec<_> = data_sets
        .iter()
        .map(|(_, data)| data)
        .collect();
    let covariance = construct_joint_covariance(&experimental_data)?;
    let inverse_covariance = invert_covariance(&covariance)?;
    let mut supports = data_sets.into_iter().map(|(support, _)| support);
    let mut support = supports.next().expect("non-empty data set list");

    for next in supports {
        if support.lx_edges != next.lx_edges || support.lq_edges != next.lq_edges {
            bail!("all jointly correlated support maps must use the same auxiliary grid");
        }
        support.weights.extend(next.weights);
    }

    let sensitivity = build_sensitivity_map(support, &inverse_covariance)?;

    normalize_sensitivity(sensitivity)
}

fn normalize_sensitivity(mut sensitivity: AuxGrid) -> Result<AuxGrid> {
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
    fn reads_cms_incjets_r04_data_and_unfolding_correlations() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("input/incjets/CMS_13TeV_R04_data_unfolding_correlated.yaml");
        let data = read_experimental_data(&path).expect("R=0.4 data should be readable");
        let covariance = construct_covariance(&data).expect("covariance should be constructible");
        invert_covariance(&covariance).expect("R=0.4 covariance should be invertible");

        assert_eq!(data.bins.len(), 78);
        assert_eq!(data.grid_bin_counts(), Some([22, 21, 19, 16].as_slice()));
        assert_eq!(data.uncertainties.len(), 29);
        assert!(data.uncertainties.iter().all(|source| source.name != "Lumi"));
        assert_eq!(covariance.shape(), (78, 78));
        assert_ne!(covariance[(0, 1)], 0.0);
    }

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
        assert_eq!(data.bins[0].value, 0.0736);
        assert_eq!(covariance.nrows(), 10);
        assert_eq!(covariance.ncols(), 10);

        let expected_first_diagonal =
            (0.00184_f64 / 0.0736).powi(2) + (0.00279_f64 / 0.0736).powi(2);
        assert!((covariance[(0, 0)] - expected_first_diagonal).abs() < 1.0e-15);
        assert_eq!(covariance[(0, 1)], 0.0);
        assert_eq!(covariance[(9, 0)], 0.0);
    }

    #[test]
    fn constructs_atlas_wplus_correlated_covariance_without_luminosity() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("input/wz/Wplus_atlas.yaml");
        let data = read_experimental_data(&path).expect("Wplus data should be readable");
        let covariance = construct_covariance(&data).expect(
            "Wplus covariance should be constructible"
        );

        assert_eq!(data.bins.len(), 11);
        assert_eq!(data.uncertainties.len(), 132);
        assert!(data.uncertainties.iter().all(|source| source.name != "AtlasLumi2011"));
        assert_eq!(covariance.nrows(), 11);
        assert_eq!(covariance.ncols(), 11);
        assert!((covariance[(0, 0)] - 3.015665465e-5).abs() < 1.0e-15);
        assert!((covariance[(0, 1)] - 2.3209e-5).abs() < 1.0e-15);
        assert!((covariance[(0, 1)] - covariance[(1, 0)]).abs() < 1.0e-15);

        invert_covariance(&covariance).expect("Wplus covariance should be invertible");
    }

    #[test]
    fn constructs_joint_atlas_wz_covariance_with_cross_process_correlations() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let wplus = read_experimental_data(&root.join("input/wz/Wplus_atlas.yaml")).unwrap();
        let wminus = read_experimental_data(&root.join("input/wz/Wminus_atlas.yaml")).unwrap();
        let zpeak = read_experimental_data(&root.join("input/wz/Zpeak_atlas.yaml")).unwrap();
        let covariance = construct_joint_covariance(&[&wplus, &wminus, &zpeak]).unwrap();

        assert_eq!((wplus.bins.len(), wminus.bins.len(), zpeak.bins.len()), (11, 11, 12));
        assert_eq!(covariance.shape(), (34, 34));
        assert!((covariance[(0, 11)] - 2.34583e-5).abs() < 1.0e-15);
        assert!((covariance[(0, 22)] - 6.1756e-6).abs() < 1.0e-15);
        assert!((covariance[(11, 22)] - 7.3797e-6).abs() < 1.0e-15);
        assert!((covariance[(0, 11)] - covariance[(11, 0)]).abs() < 1.0e-15);

        invert_covariance(&covariance).expect("joint ATLAS covariance should be invertible");
    }

    #[test]
    fn constructs_lhcb_covariances_with_uncorrelated_statistics() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));

        for (filename, number_of_bins, central, stat, sys, correlation) in [
            ("Wplus_lhcb.yaml", 8, [192.2, 178.8], [1.2, 0.9], [3.5, 3.1], -0.15),
            ("Wminus_lhcb.yaml", 8, [111.1, 104.9], [0.9, 0.7], [2.1, 1.9], -0.24),
            ("Z_lhcb.yaml", 17, [0.969, 2.84], [0.039, 0.063], [0.032, 0.05], 0.18),
        ] {
            let path = root.join("input/wz").join(filename);
            let data = read_experimental_data(&path).unwrap_or_else(|error|
                panic!("{filename} should be readable: {error}")
            );
            let covariance = construct_covariance(&data).unwrap_or_else(|error|
                panic!("{filename} covariance failed: {error}")
            );

            assert_eq!(data.bins.len(), number_of_bins);
            assert_eq!(data.uncertainties.len(), 2);
            assert_eq!(data.uncertainties[0].name, "statistical");
            assert_eq!(data.uncertainties[1].name, "systematic");
            assert_eq!(covariance.shape(), (number_of_bins, number_of_bins));

            let stat_shift: f64 = stat[0] / central[0];
            let sys_shift: f64 = sys[0] / central[0];
            let expected_diagonal = stat_shift.powi(2) + sys_shift.powi(2);
            let expected_off_diagonal =
                (((correlation * sys[0]) / central[0]) * sys[1]) / central[1];

            assert!((covariance[(0, 0)] - expected_diagonal).abs() < 1.0e-15);
            assert!((covariance[(0, 1)] - expected_off_diagonal).abs() < 1.0e-15);
            assert!((covariance[(0, 1)] - covariance[(1, 0)]).abs() < 1.0e-15);

            invert_covariance(&covariance).unwrap_or_else(|error|
                panic!("{filename} covariance is singular: {error}")
            );
        }
    }

    #[test]
    fn format_ignores_range_metadata() {
        let yaml =
            r#"
observable:
  name: "d2sigma/dxdy"
  unit: "pb"
bins:
  - ranges:
      x: { low: 0.0, high: 1.0 }
      y: { low: 2.0, high: 3.0 }
    value: 4.0
uncertainties:
  - name: "stat"
    kind: "uncorrelated"
    unit: "absolute"
    per_bin: [0.5]
"#;

        let data: ExperimentalData = serde_yaml::from_str(yaml).unwrap();
        let covariance = construct_covariance(&data).unwrap();

        assert_eq!(data.bins[0].value, 4.0);
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
            weights: vec![vec![vec![0.25, 0.75]], vec![vec![0.6, 0.4]]],
        };

        let p_matrix = build_p_matrix(&support, 2).unwrap();

        assert_eq!(p_matrix.nrows(), 2);
        assert_eq!(p_matrix.ncols(), 2);
        assert_eq!(p_matrix.column(0).into_owned(), DMatrix::from_column_slice(2, 1, &[0.25, 0.6]));
        assert_eq!(p_matrix.column(1).into_owned(), DMatrix::from_column_slice(2, 1, &[0.75, 0.4]));
    }

    #[test]
    fn builds_sensitivity_aux_grid() {
        let support = BinSupportMaps {
            lx_edges: vec![-2.0, -1.0, 0.0],
            lq_edges: vec![2.0, 3.0],
            weights: vec![vec![vec![0.25, 0.75]], vec![vec![0.6, 0.4]]],
        };
        let inverse_covariance = DMatrix::from_row_slice(2, 2, &[4.0, 0.0, 0.0, 9.0]);

        let sensitivity = build_sensitivity_map(support, &inverse_covariance).unwrap();

        assert!((sensitivity.weights[0][0] - 3.49).abs() < 1.0e-12);
        assert!((sensitivity.weights[0][1] - 3.69).abs() < 1.0e-12);
        assert_eq!(sensitivity.lx_edges, vec![-2.0, -1.0, 0.0]);
        assert_eq!(sensitivity.lq_edges, vec![2.0, 3.0]);
    }
}
