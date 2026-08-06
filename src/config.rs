//! Module to read the grids.toml file

use std::collections::HashMap;
use std::fs;
use std::path::{ Path, PathBuf };

use anyhow::{ Context, Result, anyhow };
use serde::Deserialize;

use crate::transform::GridSpec;

/// Root structure of the TOML configuration file
#[derive(Debug, Deserialize)]
pub struct GridConfig {
    path_to_grids: PathBuf,
    grids: HashMap<String, GridEntry>,
    output: OutputConfig,
    config: TransformConfig,
    selection_conf: SelectionConfig,
}

#[derive(Debug, Deserialize)]
struct TransformConfig {
    nx_x: usize,
    nx_q2: usize,
    x_min: f64,
    x_max: f64,
    q2_min: f64,
    q2_max: f64,
    threshold: f64,
}
#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum PlotStyle {
    Boundary,
    Heatmap,
}

impl Default for PlotStyle {
    fn default() -> Self {
        Self::Boundary
    }
}

#[derive(Debug, Deserialize)]
struct OutputConfig {
    path_plot: PathBuf,
    process_grid: Vec<String>,
    detail: bool,

    #[serde(default)]
    translation_ensemble: bool,

    #[serde(default)]
    groups: HashMap<String, Vec<String>>,

    #[serde(default)]
    plot_style: PlotStyle,
}

#[derive(Debug, Deserialize)]
struct SelectionConfig {
    pdgs: Option<Vec<i32>>,
    max_as: Option<u8>,
    pdf: String,
}
#[derive(Debug, Deserialize)]
struct GridEntry {
    path: PathBuf,
    label: String,
    experimental_data: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy)]
pub struct Selection {
    pdg: Option<i32>,
    max_as: Option<u8>,
}

/// Getter functions
impl GridConfig {
    pub fn plot_path(&self) -> &Path {
        &self.output.path_plot
    }

    pub fn grid_spec(&self) -> GridSpec {
        GridSpec {
            nx_x: self.config.nx_x,
            nx_q2: self.config.nx_q2,
            x_min: self.config.x_min,
            x_max: self.config.x_max,
            q2_min: self.config.q2_min,
            q2_max: self.config.q2_max,
        }
    }
    pub fn pdf(&self) -> &str {
        &self.selection_conf.pdf
    }

    pub fn threshold(&self) -> f64 {
        self.config.threshold
    }

    pub fn detail(&self) -> bool {
        self.output.detail
    }

    pub fn translation_ensemble(&self) -> bool {
        self.output.translation_ensemble
    }

    pub fn grids(&self) -> &[String] {
        &self.output.process_grid
    }

    pub fn selection(&self) -> Vec<Selection> {
        self.selection_conf.selection()
    }

    pub fn grid_label(&self, grid_name: &str) -> Result<&str> {
        self.grids
            .get(grid_name)
            .map(|entry| entry.label.as_str())
            .ok_or_else(|| anyhow!("grid '{grid_name}' is not defined in the configuration file"))
    }

    pub fn group_for_grid(&self, grid_name: &str) -> Option<&str> {
        self.output.groups.iter().find_map(|(group_name, grids)| {
            grids
                .iter()
                .any(|name| name == grid_name)
                .then_some(group_name.as_str())
        })
    }

    pub fn plot_style(&self) -> &'static str {
        match self.output.plot_style {
            PlotStyle::Boundary => "boundary",
            PlotStyle::Heatmap => "heatmap",
        }
    }

    //check if grids from [output.groups] are also in process_grid
    fn validate_groups(&self) -> Result<()> {
        let all_grouped_grids_are_processed = self.output.groups
            .values()
            .flatten()
            .all(|grid_name| self.output.process_grid.contains(grid_name));

        if !all_grouped_grids_are_processed {
            return Err(
                anyhow!("at least one grid in output.groups is not listed in output.process_grid")
            );
        }

        Ok(())
    }

    pub fn experimental_data_path(&self, grid_name: &str) -> Result<PathBuf> {
        let entry = self.grids
            .get(grid_name)
            .ok_or_else(|| {
                anyhow!("grid '{grid_name}' is not defined in the configuration file")
            })?;

        let path = entry.experimental_data
            .as_ref()
            .ok_or_else(|| { anyhow!("grid '{grid_name}' has no experimental_data entry") })?;

        Ok(path.clone())
    }
}

impl SelectionConfig {
    pub fn selection(&self) -> Vec<Selection> {
        match &self.pdgs {
            Some(pdgs) if !pdgs.is_empty() =>
                pdgs
                    .iter()
                    .copied()
                    .map(|pdg| Selection {
                        pdg: Some(pdg),
                        max_as: self.max_as,
                    })
                    .collect(),
            _ =>
                vec![Selection {
                    pdg: None,
                    max_as: self.max_as,
                }],
        }
    }
}

// if we later also want to keep always the pdg=None case (that always a reference plot is created)
// impl SelectionConfig {
//     pub fn selection(&self) -> Vec<Selection> {
//         let mut selections = vec![Selection {
//             pdg: None,
//             max_as: self.max_as,
//         }];

//         if let Some(pdgs) = &self.pdgs {
//             selections.extend(pdgs.iter().copied().map(|pdg| Selection {
//                 pdg: Some(pdg),
//                 max_as: self.max_as,
//             }));
//         }

//         selections
//     }
// }

impl Selection {
    pub fn pdg(&self) -> Option<i32> {
        self.pdg
    }

    pub fn max_as(&self) -> Option<u8> {
        self.max_as
    }
}

/// Reads and parses the TOML configuration.
pub fn read_config(config_path: &Path) -> Result<GridConfig> {
    let raw_config = fs
        ::read_to_string(config_path)
        .with_context(|| format!("could not read config file '{}'", config_path.display()))?;

    let config: GridConfig = toml
        ::from_str(&raw_config)
        .with_context(|| format!("could not parse config file '{}'", config_path.display()))?;

    config.validate_groups()?;

    Ok(config)
}

/// Resolves the full filesystem path for a configured grid name. Get the file_name for a certain grid
pub fn resolve_grid_path(config: &GridConfig, grid_name: &str) -> Result<PathBuf> {
    let file_name = config.grids
        .get(grid_name)
        .ok_or_else(|| { anyhow!("grid '{grid_name}' is not defined in the configuration file") })?;

    //check

    Ok(config.path_to_grids.join(&file_name.path))
}
// ToDo: check for correct channel masks & bins.
