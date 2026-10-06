//! Helpers to export transformed contribution data.

use std::fs::{ self, File };
use std::io::{ BufWriter, Write };
use std::path::{ Path, PathBuf };
use anyhow::{ Context, Result };
use serde::Serialize;
use serde_json::to_writer;
use tempfile::NamedTempFile;

use crate::transform::{ AuxGrid, BinSupportMaps, SummedContribution, XVar };

/// One computed distribution and the selection used to obtain it.
#[derive(Serialize)]
pub struct PlotDataset<'a> {
    pub name: &'a str,
    pub label: &'a str,
    pub group: Option<&'a str>,
    /// Names of the source grids when this distribution represents a group.
    pub members: Option<&'a [String]>,
    pub aux_grid: &'a AuxGrid,
    pub pdg: Option<i32>,
    pub leg_pdgs: Option<&'a [i32]>,
    pub channels: Option<&'a [usize]>,
    pub bins: Option<&'a [usize]>,
}

/// Shared calculation metadata and all distributions in one JSON export.
#[derive(Serialize)]
pub struct PlotExport<'a> {
    pub schema_version: u32,
    /// Calculation mode: "supp" or "sense".
    pub mode: &'a str,
    pub pdf: &'a str,
    pub pdf_member: u32,
    /// Requested order relative to the lowest alpha_s power in each source grid.
    pub max_as: Option<u8>,
    pub threshold: f64,
    /// Group weighting in support mode: "grids" or "bins".
    pub support_group_weighting: &'a str,
    pub datasets: Vec<PlotDataset<'a>>,
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

/// Write every unranked P_B,z cell using the internal logarithmic grid edges.
pub(crate) fn write_bin_support_maps_csv(
    support: &BinSupportMaps,
    file_path: &Path
) -> Result<PathBuf> {
    if let Some(parent) = file_path.parent() {
        fs
            ::create_dir_all(parent)
            .with_context(|| {
                format!("could not create output directory '{}'", parent.display())
            })?;
    }

    let file = File::create(file_path).with_context(|| {
        format!("could not create support-map file '{}'", file_path.display())
    })?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "observable_bin,iq,ix,log10_q2_low,log10_q2_high,log10_x_low,log10_x_high,p_bz"
    )?;

    for (observable_bin, weights) in support.weights.iter().enumerate() {
        for (iq, row) in weights.iter().enumerate() {
            for (ix, weight) in row.iter().enumerate() {
                writeln!(
                    writer,
                    "{},{},{},{},{},{},{},{}",
                    observable_bin,
                    iq,
                    ix,
                    support.lq_edges[iq],
                    support.lq_edges[iq + 1],
                    support.lx_edges[ix],
                    support.lx_edges[ix + 1],
                    weight
                )?;
            }
        }
    }

    writer
        .flush()
        .with_context(|| {
            format!("could not flush support-map file '{}'", file_path.display())
        })?;

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

/// Function to create tmp file where the auxilary grid data is written into a json file.
/// The json file can then later be processed by the plotting script
pub fn write_auxgrids_tmp_json(auxgrid: &[PlotDataset<'_>]) -> Result<NamedTempFile> {
    let mut file = NamedTempFile::new().context(
        "could not create temporary file for auxiliary grid JSON"
    )?;

    to_writer(file.as_file_mut(), auxgrid).context(
        "could not serialize auxiliary grid to temporary JSON"
    )?;
    file.as_file_mut().flush().context("could not flush temporary auxiliary grid JSON")?;

    Ok(file)
}
