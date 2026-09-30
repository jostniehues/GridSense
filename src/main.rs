use std::{ collections::BTreeMap, fs::File };
use std::path::PathBuf;

use anyhow::{ Context, Result };
use clap::{ Parser, ValueEnum };
use pineappl::grid::Grid;

mod config;
mod grids;
mod transform;
mod export;
mod sensitivity;
use config::{ read_config, resolve_grid_paths, SupportGroupWeighting };
use std::io::{ self, Write };
use grids::extract_kinematics;
use std::path::Path;
use std::process::Command;

use export::{
    write_auxgrids_tmp_json,
    write_bin_support_maps_csv,
    write_q2xcontrib_csv,
    PlotDataset,
};
use transform::{
    average_support_maps,
    build_per_bin_support_maps,
    ranked_cells,
    sum_over_x,
    BinSupportMaps,
};
use sensitivity::{
    calculate_joint_sensitivity_map,
    calculate_sensitivity_map,
    read_experimental_data,
    ExperimentalData,
};

use crate::{ config::{ GridConfig, Selection }, transform::AuxGrid };

struct PlotGrid {
    name: String,
    label: String,
    group: Option<String>,
    pdg: Option<i32>,
    aux_grid: AuxGrid,
    max_as: Option<u8>,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Mode {
    Supp,
    Sense,
}

#[derive(Debug, Parser)]
#[command(
    name = "gridsense",
    version,
    about = "Build kinematic support from a PineAPPL interpolation grid"
)]
struct Cli {
    /// Path to the TOML file that contains the grid directory and mappings.
    #[arg(short, long, default_value = "grids.toml")]
    config: PathBuf,

    /// Operation mode: kinematic support or uncertainty-weighted sensitivity.
    #[arg(long, value_enum, default_value = "supp")]
    mode: Mode,
}

/// Load one grid and build its per-observable-bin support maps.
/// In the case of mutlitple grids, concatenate their observable bins
fn load_bin_support_grid(
    grid_name: &str,
    config: &GridConfig,
    selection: &Selection,
    write_contributions: bool
) -> Result<BinSupportMaps> {
    load_bin_support_grid_with_counts(grid_name, config, selection, write_contributions, None)
}

fn load_bin_support_grid_with_counts(
    grid_name: &str,
    config: &GridConfig,
    selection: &Selection,
    write_contributions: bool,
    grid_bin_counts: Option<&[usize]>
) -> Result<BinSupportMaps> {
    let grid_paths = resolve_grid_paths(config, grid_name)?;

    if let Some(counts) = grid_bin_counts {
        if counts.len() != grid_paths.len() {
            anyhow::bail!(
                "experimental data for '{grid_name}' defines {} grid_bin_counts entries for {} grid paths",
                counts.len(),
                grid_paths.len()
            );
        }
    }

    let mut combined_kinematics = Vec::new();

    for (grid_index, grid_path) in grid_paths.into_iter().enumerate() {
        let grid_file = File::open(&grid_path).with_context(||
            format!("could not open grid '{}'", grid_path.display())
        )?;
        let grid = Grid::read(grid_file).with_context(||
            format!("could not read grid '{}'", grid_path.display())
        )?;

        let mut kinematics = extract_kinematics(&grid, selection, config.pdf()).with_context(|| {
            format!("could not extract kinematics from grid '{}'", grid_path.display())
        })?;

        if let Some(number_of_bins) = grid_bin_counts.map(|counts| counts[grid_index]) {
            if number_of_bins > kinematics.len() {
                anyhow::bail!(
                    "experimental data requests {number_of_bins} bins from grid '{}', but it contains only {}",
                    grid_path.display(),
                    kinematics.len()
                );
            }

            kinematics.truncate(number_of_bins);
        }

        combined_kinematics.extend(kinematics);
    }

    if config.detail() && write_contributions {
        let summed_points_x1 = sum_over_x(&combined_kinematics, transform::XVar::X1)?;
        let summed_points_x2 = sum_over_x(&combined_kinematics, transform::XVar::X2)?;
        let output_dir = Path::new("output/data/contrib");
        let path = output_dir.join(format!("{grid_name}.csv"));

        write_q2xcontrib_csv(&summed_points_x1, &summed_points_x2, &path)?;
    }

    if let Some(bins) = selection.bins() {
        let number_bins = combined_kinematics.len();

        for &bin in bins {
            if bin >= number_bins {
                anyhow::bail!(
                    "bin index {bin} is out of range for '{grid_name}': \
                 the grid contains {number_bins} bins"
                );
            }
        }
        combined_kinematics = combined_kinematics
            .into_iter()
            .enumerate()
            .filter(|(index, _)| bins.binary_search(index).is_ok())
            .map(|(_, contributions)| contributions)
            .collect();
    }

    // Build the normalized support map from the three-dimensional auxiliary binning.
    let grid_spec = config.grid_spec();

    let per_bin_support = build_per_bin_support_maps(&combined_kinematics, &grid_spec)?;

    if config.detail() {
        let path = Path::new("output/data/support").join(format!("{grid_name}.csv"));
        write_bin_support_maps_csv(&per_bin_support, &path)?;
    }
    Ok(per_bin_support)
}

fn combine_grouped_support(
    support_maps: Vec<BinSupportMaps>,
    weighting: SupportGroupWeighting
) -> Result<AuxGrid> {
    let Some(reference) = support_maps.first() else {
        anyhow::bail!("cannot combine an empty support group");
    };

    if
        support_maps
            .iter()
            .any(|support| {
                support.lx_edges != reference.lx_edges || support.lq_edges != reference.lq_edges
            })
    {
        anyhow::bail!("cannot combine support maps with different auxiliary grids");
    }

    let number_of_grids = support_maps.len();
    let mut support_maps = support_maps.into_iter();
    let first = support_maps.next().expect("support group was checked as non-empty");

    match weighting {
        SupportGroupWeighting::Grids => {
            let mut combined = average_support_maps(first);

            for support in support_maps {
                combined.merge(average_support_maps(support));
            }

            for weight in combined.weights.iter_mut().flatten() {
                *weight /= number_of_grids as f64;
            }

            Ok(combined)
        }
        SupportGroupWeighting::Bins => {
            let mut combined = first;

            for support in support_maps {
                combined.weights.extend(support.weights);
            }

            Ok(average_support_maps(combined))
        }
    }
}

/// Let user confirm that (s)he wants to override the path if it already exsits
fn confirm_overwrite(path: &Path) -> Result<bool> {
    if path.exists() {
        print!("File '{}' already exists. Overwrite? [y/n]: ", path.display());
        io::stdout().flush()?;

        let mut input = String::new();
        io::stdin().read_line(&mut input)?;

        let input = input.trim().to_lowercase();

        Ok(input == "y" || input == "yes" || input == "j" || input == "ja")
    } else {
        Ok(true)
    }
}

// Load one grid and its experimental input for a sensitivity calculation.
fn load_sensitivity_grid(
    grid_name: &str,
    experimental_path: &Path,
    config: &GridConfig,
    selection: &Selection
) -> Result<(BinSupportMaps, ExperimentalData)> {
    let experimental_data = read_experimental_data(experimental_path)?;
    let support = load_bin_support_grid_with_counts(
        grid_name,
        config,
        selection,
        false,
        experimental_data.grid_bin_counts()
    )?;

    Ok((support, experimental_data))
}

fn main() {
    if let Err(error) = run() {
        eprintln!("Error: {error}");
        std::process::exit(1);
    }
}

/// Executes the CLI workflow and returns a readable error on failure.
fn run() -> Result<()> {
    let cli = Cli::parse();
    let config = read_config(&cli.config)?;

    match cli.mode {
        Mode::Supp => run_support(&config),
        Mode::Sense => run_sensitivity(&config),
    }
}

/// Export the calculated auxiliary grids and pass them to the Python plotter.
fn plot_aux_grids(config: &GridConfig, aux_grids: &[PlotGrid], mode: Mode) -> Result<()> {
    let plot_path = config.plot_path();
    let plot_datasets: Vec<_> = aux_grids
        .iter()
        .map(|grid| PlotDataset {
            name: &grid.name,
            label: &grid.label,
            group: grid.group.as_deref(),
            pdg: &grid.pdg,
            aux_grid: &grid.aux_grid,
            max_as: &grid.max_as,
        })
        .collect();

    let tmp_file = write_auxgrids_tmp_json(&plot_datasets)?;

    if let Some(parent) = plot_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    let status = Command::new("python3")
        .arg("python/plot_auxgrid.py")
        .arg(tmp_file.path())
        .arg(plot_path)
        .arg(config.plot_style())
        .arg(match mode {
            Mode::Supp => "supp",
            Mode::Sense => "sense",
        })
        .arg(config.subtitle().unwrap_or(""))
        .status()
        .context("could not start python plotting script")?;

    if !status.success() {
        anyhow::bail!("python plotting script failed with status {status}");
    }

    Ok(())
}

fn run_support(config: &GridConfig) -> Result<()> {
    // final plot path
    let plot_path = config.plot_path();

    let overwrite = confirm_overwrite(plot_path)?;
    //if user does not want to overwrite
    if !overwrite {
        println!("Aborted, existing plot was not overwritten.");
        return Ok(());
    }

    // create temporary json file
    let mut aux_grids = Vec::new();

    for selec in config.selection() {
        let mut grouped_support = BTreeMap::<String, Vec<BinSupportMaps>>::new();

        for grid_name in config.grids() {
            let bin_support = load_bin_support_grid(grid_name, config, &selec, true)?;
            let label = config.grid_label(grid_name)?;

            if let Some(group_name) = config.group_for_grid(grid_name) {
                grouped_support.entry(group_name.to_string()).or_default().push(bin_support);
            } else {
                let aux_grid = average_support_maps(bin_support);
                aux_grids.push(PlotGrid {
                    name: grid_name.to_string(),
                    label: label.to_string(),
                    group: None,
                    pdg: selec.pdg(),
                    aux_grid: ranked_cells(aux_grid, config.threshold())?,
                    max_as: selec.max_as(),
                });
            }
        }

        for (group_name, support_maps) in grouped_support {
            let aux_grid = combine_grouped_support(support_maps, config.support_group_weighting())?;

            aux_grids.push(PlotGrid {
                name: group_name.clone(),
                label: group_name.clone(),
                group: Some(group_name),
                pdg: selec.pdg(),
                aux_grid: ranked_cells(aux_grid, config.threshold())?,
                max_as: selec.max_as(),
            });
        }
    }

    plot_aux_grids(config, &aux_grids, Mode::Supp)
}

//Todo: avoid code duplication
fn run_sensitivity(config: &GridConfig) -> Result<()> {
    //bin selection not implemented for sensititivy mode
    if config.selection().iter().any(|selection| selection.bins().is_some()) {
        anyhow::bail!(
            "bin selection is currently supported only in supp mode"
        );
    }

    let plot_path = config.plot_path();

    let overwrite = confirm_overwrite(plot_path)?;
    //if user does not want to overwrite
    if !overwrite {
        println!("Aborted, existing plot was not overwritten.");
        return Ok(());
    }

    // create temporary json file
    let mut aux_grids = Vec::new();

    for selection in config.selection() {
        let mut grouped_data = BTreeMap::<String, Vec<(BinSupportMaps, ExperimentalData)>>::new();

        for grid_name in config.grids() {
            let experimental_path = config.experimental_data_path(grid_name)?;
            let data_set = load_sensitivity_grid(
                grid_name,
                &experimental_path,
                config,
                &selection
            )?;

            if let Some(group_name) = config.group_for_grid(grid_name) {
                grouped_data.entry(group_name.to_string()).or_default().push(data_set);
            } else {
                let (support, experimental_data) = data_set;
                let sensitivity = calculate_sensitivity_map(support, &experimental_data)?;
                let aux_grid = ranked_cells(sensitivity, config.threshold())?;
                let label = config.grid_label(grid_name)?;

                aux_grids.push(PlotGrid {
                    name: grid_name.to_string(),
                    label: label.to_string(),
                    group: None,
                    pdg: selection.pdg(),
                    aux_grid,
                    max_as: selection.max_as(),
                });
            }
        }

        for (group_name, data_sets) in grouped_data {
            let sensitivity = calculate_joint_sensitivity_map(data_sets)?;
            let aux_grid = ranked_cells(sensitivity, config.threshold())?;

            aux_grids.push(PlotGrid {
                name: group_name.clone(),
                label: group_name.clone(),
                group: Some(group_name),
                pdg: selection.pdg(),
                aux_grid,
                max_as: selection.max_as(),
            });
        }
    }

    plot_aux_grids(config, &aux_grids, Mode::Sense)
}
