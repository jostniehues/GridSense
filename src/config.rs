//! Module to read the grids.toml file

use std::collections::{ BTreeSet, HashMap };
use std::fs;
use std::path::{ Path, PathBuf };

use anyhow::{ Context, Result, anyhow };
use serde::{ Deserialize, Deserializer, de };

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
    groups: HashMap<String, Vec<String>>,

    #[serde(default)]
    plot_style: PlotStyle,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct SelectionConfig {
    pdgs: Option<Vec<i32>>,
    #[serde(default)]
    leg_pdgs: Option<Vec<i32>>,
    max_as: Option<u8>,
    #[serde(default, deserialize_with = "deserialize_channels")]
    channels: Option<Vec<usize>>,
    pdf: String,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum ChannelItem {
    Index(usize),
    Range(String),
}

fn deserialize_channels<'de, D>(deserializer: D) -> Result<Option<Vec<usize>>, D::Error>
    where D: Deserializer<'de>
{
    let Some(items) = Option::<Vec<ChannelItem>>::deserialize(deserializer)? else {
        return Ok(None);
    };

    if items.is_empty() {
        return Err(
            de::Error::custom("channels must not be empty; omit the field to select all channels")
        );
    }

    let mut channels = BTreeSet::new();

    for item in items {
        match item {
            ChannelItem::Index(index) => {
                channels.insert(index);
            }
            ChannelItem::Range(value) => {
                let value = value.trim();

                if let Some((start, end)) = value.split_once('-') {
                    let start = start
                        .trim()
                        .parse::<usize>()
                        .map_err(|_| {
                            de::Error::custom(format!("invalid channel range '{value}'"))
                        })?;
                    let end = end
                        .trim()
                        .parse::<usize>()
                        .map_err(|_| {
                            de::Error::custom(format!("invalid channel range '{value}'"))
                        })?;

                    if start > end {
                        return Err(
                            de::Error::custom(
                                format!("channel range '{value}' has decreasing bounds")
                            )
                        );
                    }

                    channels.extend(start..=end);
                } else {
                    let index = value
                        .parse::<usize>()
                        .map_err(|_| {
                            de::Error::custom(format!("invalid channel index '{value}'"))
                        })?;
                    channels.insert(index);
                }
            }
        }
    }

    Ok(Some(channels.into_iter().collect()))
}
#[derive(Debug, Deserialize)]
struct GridEntry {
    path: PathBuf,
    label: String,
    experimental_data: Option<PathBuf>,
}

#[derive(Debug, Clone)]
pub struct Selection {
    pdg: Option<i32>,
    max_as: Option<u8>,
    leg_pdgs: Option<Vec<i32>>,
    channels: Option<Vec<usize>>,
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
                        leg_pdgs: None,
                        max_as: self.max_as,
                        channels: None,
                    })
                    .collect(),

            _ =>
                vec![Selection {
                    pdg: None,
                    leg_pdgs: self.leg_pdgs.clone(),
                    max_as: self.max_as,
                    channels: self.channels.clone(),
                }],
        }
    }

    fn validate(&self) -> Result<()> {
        if self.pdgs.as_ref().is_some_and(Vec::is_empty) {
            anyhow::bail!("selection_conf.pdgs must not be empty; omit the field instead");
        }

        if self.leg_pdgs.as_ref().is_some_and(Vec::is_empty) {
            anyhow::bail!("selection_conf.leg_pdgs must not be empty; omit the field instead");
        }

        if self.pdgs.is_some() && self.channels.is_some() {
            anyhow::bail!(
                "selection_conf.pdgs and selection_conf.channels cannot be used together"
            );
        }

        if self.pdgs.is_some() && self.leg_pdgs.is_some() {
            anyhow::bail!(
                "selection_conf.pdgs and selection_conf.leg_pdgs cannot be used together"
            );
        }

        Ok(())
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

    pub fn channels(&self) -> Option<&[usize]> {
        self.channels.as_deref()
    }
    pub fn leg_pdgs(&self) -> Option<&[i32]> {
        self.leg_pdgs.as_deref()
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

    // add central method later when there are more validation steps
    config.selection_conf.validate()?;
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_channel_indices_and_ranges() {
        let config: SelectionConfig = toml
            ::from_str(r#"
pdf = "CT18NNLO"
channels = ["1-5", 8, 11]
"#)
            .unwrap();

        assert_eq!(config.channels, Some(vec![1, 2, 3, 4, 5, 8, 11]));
    }

    #[test]
    fn rejects_decreasing_channel_ranges() {
        let error = toml
            ::from_str::<SelectionConfig>(r#"
pdf = "CT18NNLO"
channels = ["5-1"]
"#)
            .unwrap_err();

        assert!(error.to_string().contains("decreasing bounds"));
    }
}
