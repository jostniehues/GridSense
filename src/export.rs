//! Helpers to export transformed contribution data.

use std::fs::{ self, File };
use std::io::{ BufWriter, Write };
use std::path::{ Path, PathBuf };
use anyhow::{ Context, Result };
use serde::Serialize;
use serde_json::to_writer;
use tempfile::NamedTempFile;

use crate::transform::{ AuxGrid, SummedContribution, XVar };

#[derive(Serialize)]
pub struct PlotDataset<'a> {
    pub name: &'a str,
    pub label: &'a str,
    pub group: Option<&'a str>,
    pub aux_grid: &'a AuxGrid,
    pub pdg: &'a Option<i32>,
    pub max_as: &'a Option<u8>,
}

/// Writes (bin, q2, x, xvar, contrib) rows to a CSV file.
/// Therefore it needs the summed contribution over x1 respective x2
pub fn write_q2xcontrib_csv(
    points_x1: &[Vec<SummedContribution>],
    points_x2: &[Vec<SummedContribution>],
    file_path: &Path
) -> Result<PathBuf> {
    if let Some(parent) = file_path.parent() {
        fs
            ::create_dir_all(parent)
            .with_context(|| {
                format!("could not create output directory '{}'", parent.display())
            })?;
    }

    let file = File::create(file_path).with_context(||
        format!("could not create file '{}'", file_path.display())
    )?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "bin,q2,x,xvar,contrib").with_context(||
        format!("could not write header to '{}'", file_path.display())
    )?;

    write_projection_rows(&mut writer, points_x1, file_path)?;
    write_projection_rows(&mut writer, points_x2, file_path)?;

    writer
        .flush()
        .with_context(|| format!("could not flush writer for '{}'", file_path.display()))?;

    Ok(file_path.to_path_buf())
}

fn xvar_name(xvar: XVar) -> &'static str {
    match xvar {
        XVar::X1 => "x1",
        XVar::X2 => "x2",
    }
}

fn write_projection_rows(
    writer: &mut BufWriter<File>,
    points: &[Vec<SummedContribution>],
    file_path: &Path
) -> Result<()> {
    for (bin_index, bin_points) in points.iter().enumerate() {
        for point in bin_points {
            writeln!(
                writer,
                "{},{},{},{},{}",
                bin_index,
                point.q2,
                point.x,
                xvar_name(point.xvar),
                point.value
            ).with_context(|| format!("could not write row to '{}'", file_path.display()))?;
        }
    }

    Ok(())
}

/// Function in order to create tmp file where the auxilary grid data is written into a json file.
/// The json file can then later be processed by the plotting script
pub fn write_auxgrids_tmp_json(auxgrid: &[PlotDataset<'_>]) -> Result<NamedTempFile> {
    let mut file = NamedTempFile::new()
        .context("could not create temporary file for auxiliary grid JSON")?;

    to_writer(file.as_file_mut(), auxgrid)
        .context("could not serialize auxiliary grid to temporary JSON")?;
    file
        .as_file_mut()
        .flush()
        .context("could not flush temporary auxiliary grid JSON")?;

    Ok(file)
}
