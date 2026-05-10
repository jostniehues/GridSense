use std::{ fs::File };
use std::path::PathBuf;

use anyhow::{ Context, Result };
use clap::Parser;
use pineappl::grid::Grid;

mod config;
mod grids;
mod transform;
mod export;
use config::{ read_config, resolve_grid_path };
use std::io::{ self, Write };
use grids::extract_kinematics;
use std::path::Path;
use std::process::Command;

use export::{ write_auxgrids_tmp_json, write_q2xcontrib_csv, PlotDataset };
use transform::{ build_relevance, ranked_cells, sort_into_auxgrid, sum_over_x,union_aux_grid };

use crate::{ config::{ GridConfig, Selection }, transform::AuxGrid };

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
}

/// helper function to process multiple grids
fn process_grid(grid_name: &str, config: &GridConfig, selection: &Selection) -> Result<AuxGrid> {
    let grid_path = resolve_grid_path(&config, grid_name)?;

    let grid_file = File::open(&grid_path).with_context(||
        format!("could not open grid '{}'", grid_path.display())
    )?;
    let grid = Grid::read(grid_file).with_context(||
        format!("could not read grid '{}'", grid_path.display())
    )?;

    let kinematics = extract_kinematics(&grid, selection)?;
    let summed_points_x1 = sum_over_x(&kinematics, transform::XVar::X1)?;
    let summed_points_x2 = sum_over_x(&kinematics, transform::XVar::X2)?;

    //TODO: Dieses Spezifkation soll geprinted werden, wenn detail flag in toml gegeben ist. Derzeit aber hardgecoded
    let path = Path::new("output/data/grid.csv");
    write_q2xcontrib_csv(&summed_points_x1, &summed_points_x2, path)?;

    // Grid operations
    let relevance_x1 = build_relevance(&summed_points_x1)?;
    let relevance_x2 = build_relevance(&summed_points_x2)?;
    let grid_spec = config.grid_spec();
    let aux_grid_x1 = sort_into_auxgrid(&relevance_x1, &grid_spec)?;
    let aux_grid_x2 = sort_into_auxgrid(&relevance_x2, &grid_spec)?;
    let ranked_aux_grid_x1 = ranked_cells(aux_grid_x1, config.threshold())?;
    let ranked_aux_grid_x2 = ranked_cells(aux_grid_x2, config.threshold())?;
    let union_grid = union_aux_grid(&ranked_aux_grid_x1, &ranked_aux_grid_x2)?;

    Ok(union_grid)
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
        for grid_name in config.grids() {
            let aux_grid = process_grid(grid_name, &config, &selec)?;
            let label = config.grid_label(grid_name)?;
            aux_grids.push((grid_name.to_string(), label.to_string(), selec.pdg(), aux_grid));
        }
    }

    let plot_datasets: Vec<_> = aux_grids
        .iter()
        .map(|(name, label, pdg, aux_grid)| PlotDataset {
            name,
            label,
            pdg,
            aux_grid,
        })
        .collect();

    let tmp_file = write_auxgrids_tmp_json(&plot_datasets)?;


    if let Some(parent) = plot_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    // trigger python
    let status = Command::new("python3")
        .arg("python/plot_auxgrid.py")
        .arg(tmp_file.path())
        .arg(plot_path)
        .status()
        .with_context(|| "could not start python plotting script")?;

    if !status.success() {
        anyhow::bail!("python plotting script failed with status {status}");
    }

    Ok(())
}
