//! Module to extract relevant data from pineappl grids
//!
use anyhow::{ Context, Ok, Result };
use lhapdf::Pdf;
use pineappl::{
    convolutions::{ Conv, ConvType, ConvolutionCache },
    grid::{ Grid, KinematicsContribution },
};
use crate::config::Selection;
use pineappl::evolution::gluon_has_pid_zero;

pub fn extract_kinematics(
    grid: &Grid,
    selection: &Selection
) -> Result<Vec<Vec<KinematicsContribution>>> {
    lhapdf::set_verbosity(0);

    let pdf = Pdf::with_setname_and_member("CT18NNLO", 0).context(
        "could not load LHAPDF set 'CT18NNLO' member 0"
    )?;
    let pdf_alphas = Pdf::with_setname_and_member("CT18NNLO", 0).context(
        "could not load LHAPDF set 'CT18NNLO' member 0 for alpha_s"
    )?;
    let mut xfx = |id, x, q2| pdf.xfx_q2(id, x, q2);
    let mut alphas = |q2: f64| pdf_alphas.alphas_q2(q2);

    let mut convolution_cache = ConvolutionCache::new(
        vec![Conv::new(ConvType::UnpolPDF, 2212)],
        vec![&mut xfx],
        &mut alphas
    );

    let channel_mask = vec![true; 116];
    //channel_mask[0] = true;

    let order_mask = build_order_mask(grid, selection);
    dbg!(&order_mask);

 
    //some grids may use 0 or 21 for gluon pdg
    let pdg = match selection.pdg() {
        Some(0) if gluon_has_pid_zero(grid) => Some(0),
        Some(0) => Some(21),
        other => other,
    };

  
    // "regular" convolution as a sanity check
    // TBD
    let _regular_conv = grid.convolve(
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
        &[],
        &[(1.0, 1.0, 1.0)],
        pdg
    );

    Ok(contributions)
}


// Build a order mask matching the condition in selection since different bins have different available order. 
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

