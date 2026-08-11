//! Module to extract relevant data from pineappl grids
//!
use anyhow::{ Context, Ok, Result };
use lhapdf::Pdf;
use pineappl::{
    convolutions::{ Conv, ConvType, ConvolutionCache },
    grid::{ Grid, KinematicsContribution, PdgFilter },
};
use crate::config::Selection;
use pineappl::evolution::gluon_has_pid_zero;

//some grids may use 0 or 21 for gluon pdg
fn translate_pdg(grid: &Grid, pdg: i32) -> i32 {
    match pdg {
        0 if gluon_has_pid_zero(grid) => 0,
        0 => 21,
        other => other,
    }
}

pub fn extract_kinematics(
    grid: &Grid,
    selection: &Selection,
    pdf_name: &str
) -> Result<Vec<Vec<KinematicsContribution>>> {
    lhapdf::set_verbosity(0);
    let pdf = Pdf::with_setname_and_member(pdf_name, 0).with_context(|| {
        format!("could not load LHAPDF set '{pdf_name}' member 0")
    })?;

    let pdf_alphas = Pdf::with_setname_and_member(pdf_name, 0).with_context(|| {
        format!("could not load LHAPDF set '{pdf_name}' member 0 for alpha_s")
    })?;

    let mut xfx = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut alphas = |q2: f64| pdf_alphas.alphas_q2(q2);

    let mut convolution_cache = ConvolutionCache::new(
        vec![Conv::new(ConvType::UnpolPDF, 2212)],
        vec![&mut xfx],
        &mut alphas
    );

    let order_mask = build_order_mask(grid, selection);
    let channel_mask = build_channel_mask(grid, selection)?;

    let pdg_filter = if let Some(pdg) = selection.pdg() {
        PdgFilter::InitialState(translate_pdg(grid, pdg))
    } else if let Some(pdgs) = selection.leg_pdgs() {
        PdgFilter::MatchingLegs(
            pdgs
                .iter()
                .copied()
                .map(|pdg| translate_pdg(grid, pdg))
                .collect()
        )
    } else {
        PdgFilter::None
    };

    let regular_conv = grid.convolve(
        &mut convolution_cache,
        &order_mask,
        &[],
        &channel_mask,
        &[(1.0, 1.0, 1.0)]
    );

    let contributions = grid.convolve_extract_kinematics(
        &mut convolution_cache,
        &order_mask,
        &[],
        &channel_mask,
        &[(1.0, 1.0, 1.0)],
        &pdg_filter
    );

    //closure test
    if matches!(pdg_filter, PdgFilter::None) {
        check_convolution_closure(&regular_conv, &contributions)?;
    }

    Ok(contributions)
}

/// Build a mask for the zero-based PineAPPL channel indices selected in the
/// configuration. An empty mask means that PineAPPL uses all channels.
fn build_channel_mask(grid: &Grid, selection: &Selection) -> Result<Vec<bool>> {
    let Some(channels) = selection.channels() else {
        return Ok(Vec::new());
    };

    let number_of_channels = grid.channels().len();
    let mut mask = vec![false; number_of_channels];

    for &channel in channels {
        if channel >= number_of_channels {
            anyhow::bail!(
                "channel index {channel} is out of range for a grid with {number_of_channels} channels (valid indices: 0..{})",
                number_of_channels.saturating_sub(1)
            );
        }

        mask[channel] = true;
    }

    Ok(mask)
}

/// Build a order mask matching the condition in selection since different bins have different available order.
fn build_order_mask(grid: &Grid, selection: &Selection) -> Vec<bool> {
    grid.orders()
        .iter()
        .map(|order| {
            let as_ok = selection.max_as().is_none_or(|max_as| order.alphas <= max_as);
            let no_scale_logs = order.logxif == 0 && order.logxir == 0 && order.logxia == 0;

            as_ok && no_scale_logs
        })
        .collect()
}

/// Function to check that the convolution over every single grid point reproduces the pineappl convolution
fn check_convolution_closure(
    regular: &[f64],
    node_contrib: &[Vec<KinematicsContribution>]
) -> Result<()> {
    if regular.len() != node_contrib.len() {
        anyhow::bail!(
            "sanity check failed: convolution and extracted contributions have different numbers of bins"
        );
    }
    const REL_TOLERANCE: f64 = 1.0e-8;
    const ABS_TOLERANCE: f64 = 1.0e-12;

    for (bin_index, (expected, bin)) in regular.iter().zip(node_contrib).enumerate() {
        let reconstructed: f64 = bin
            .iter()
            .map(|point| point.value)
            .sum(); //sum over all interpolation grid points

        let difference = (reconstructed - expected).abs();
        let tolerance = ABS_TOLERANCE + REL_TOLERANCE * expected.abs();

        if difference > tolerance {
            anyhow::bail!(
                "sanity check failed in observable bin {bin_index}: \
                 regular convolution = {expected}, \
                 reconstructed contribution = {reconstructed}"
            );
        }
    }
    Ok(())
}
