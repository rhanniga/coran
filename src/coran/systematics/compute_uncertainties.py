"""
Compute systematic uncertainties from variation analysis results.

This script reads the analysis_out.root file containing central values and variations,
and computes systematic uncertainties using RMS method.
"""

import ROOT as rt
import numpy as np
from pathlib import Path


def compute_rms_uncertainty(values: list[float]) -> float:
    """
    Compute RMS of a list of values.
    
    Args:
        values: List of relative differences abs(1 - variation/central)
    
    Returns:
        RMS value
    """
    if not values:
        return 0.0
    return np.sqrt(np.mean(np.array(values) ** 2))


def get_variation_type(variation_name: str) -> str | None:
    """
    Extract the variation type from the variation name.
    
    Args:
        variation_name: Full variation name (e.g., "variation_signal_wide")
    
    Returns:
        Variation type (e.g., "signal", "sideband") or None if central
    """
    if variation_name == "central":
        return None
    
    # Expected format: "variation_{type}_{name}"
    parts = variation_name.split("_")
    if len(parts) >= 2 and parts[0] == "variation":
        return parts[1]
    
    return None


def compute_deltaphi_systematics(
    central_hist: rt.TH1D,
    variation_hists: dict[str, dict[str, rt.TH1D]],
) -> tuple[rt.TH1D, rt.TH1D, dict[str, float]]:
    """
    Compute systematic uncertainties for delta phi distribution.
    
    For each variation type, computes RMS of abs(1 - variation/central) across
    all delta phi bins, then combines in quadrature.
    
    Args:
        central_hist: Central delta phi histogram
        variation_hists: Dict of {variation_type: {variation_name: histogram}}
    
    Returns:
        Tuple of (hist_with_stat_errors, hist_with_sys_errors, sys_by_type)
    """
    n_bins = central_hist.GetNbinsX()
    
    # Create output histograms
    hist_stat = central_hist.Clone("delta_phi_hist_stat")
    hist_sys = central_hist.Clone("delta_phi_hist_sys")
    
    # Track systematics by variation type for bookkeeping
    systematics_by_type = {}
    
    # Compute systematic uncertainty for each bin
    for bin_i in range(1, n_bins + 1):
        central_value = central_hist.GetBinContent(bin_i)
        
        # Skip bins with zero or negative central values
        if central_value <= 0:
            hist_sys.SetBinError(bin_i, 0)
            continue
        
        # Collect relative differences by variation type
        rel_diffs_by_type = {}
        
        for var_type, var_hists in variation_hists.items():
            rel_diffs = []
            
            for var_name, var_hist in var_hists.items():
                var_value = var_hist.GetBinContent(bin_i)
                
                if var_value > 0:
                    rel_diff = abs(1.0 - var_value / central_value)
                    rel_diffs.append(rel_diff)
            
            if rel_diffs:
                rms = compute_rms_uncertainty(rel_diffs)
                rel_diffs_by_type[var_type] = rms
        
        # Combine systematics from different variation types in quadrature
        total_sys_rel = np.sqrt(sum(rms**2 for rms in rel_diffs_by_type.values()))
        total_sys_abs = total_sys_rel * central_value
        
        hist_sys.SetBinError(bin_i, total_sys_abs)
    
    # Compute bin-averaged systematic uncertainty by type for bookkeeping
    for var_type in variation_hists.keys():
        bin_sys = []
        for bin_i in range(1, n_bins + 1):
            central_value = central_hist.GetBinContent(bin_i)
            if central_value <= 0:
                continue
            
            rel_diffs = []
            for var_name, var_hist in variation_hists[var_type].items():
                var_value = var_hist.GetBinContent(bin_i)
                if var_value > 0:
                    rel_diff = abs(1.0 - var_value / central_value)
                    rel_diffs.append(rel_diff)
            
            if rel_diffs:
                bin_sys.append(compute_rms_uncertainty(rel_diffs))
        
        if bin_sys:
            systematics_by_type[var_type] = np.mean(bin_sys)
    
    return hist_stat, hist_sys, systematics_by_type


def compute_scalar_systematics(
    central_value: float,
    variation_values: dict[str, dict[str, float]],
) -> tuple[float, dict[str, float]]:
    """
    Compute systematic uncertainty for a scalar quantity (yield or width).
    
    For each variation type, computes RMS of abs(1 - variation/central),
    then combines in quadrature.
    
    Args:
        central_value: Central value of the quantity
        variation_values: Dict of {variation_type: {variation_name: value}}
    
    Returns:
        Tuple of (total_systematic_uncertainty, sys_by_type)
    """
    if central_value <= 0:
        return 0.0, {}
    
    systematics_by_type = {}
    
    for var_type, var_values in variation_values.items():
        rel_diffs = []
        
        for var_name, var_value in var_values.items():
            if var_value > 0:
                rel_diff = abs(1.0 - var_value / central_value)
                rel_diffs.append(rel_diff)
        
        if rel_diffs:
            rms = compute_rms_uncertainty(rel_diffs)
            systematics_by_type[var_type] = rms
    
    # Combine in quadrature
    total_sys_rel = np.sqrt(sum(rms**2 for rms in systematics_by_type.values()))
    total_sys_abs = total_sys_rel * central_value
    
    return total_sys_abs, systematics_by_type


def process_mult_pt_bin(
    input_file: rt.TFile,
    output_file: rt.TFile,
    mult_bin_name: str,
    pt_bin_name: str,
):
    """
    Process one multiplicity/pt bin and compute systematic uncertainties.
    
    Args:
        input_file: Input ROOT file with analysis results
        output_file: Output ROOT file for results with systematics
        mult_bin_name: Name of multiplicity bin
        pt_bin_name: Name of pt bin
    """
    # Get central values
    central_dir_path = f"{mult_bin_name}/{pt_bin_name}/central"
    central_dir = input_file.GetDirectory(central_dir_path)
    
    if not central_dir:
        print(f"Warning: Central directory not found for {mult_bin_name}/{pt_bin_name}")
        return
    
    # Read central histogram and quantities
    central_hist = central_dir.Get("delta_phi_hist")
    if not central_hist:
        print(f"Warning: Central histogram not found for {mult_bin_name}/{pt_bin_name}")
        return
    
    central_hist = central_hist.Clone()
    central_hist.SetDirectory(0)
    
    # Read central scalar quantities
    central_quantities = {}
    quantity_names = [
        "ns_width", "as_width",
        "ns_yield", "as_yield", "ue_yield", "total_yield"
    ]
    
    for qty_name in quantity_names:
        qty_param = central_dir.Get(qty_name)
        if qty_param:
            central_quantities[qty_name] = qty_param.GetVal()
        
        # Also get statistical errors
        qty_error_param = central_dir.Get(f"{qty_name}_error")
        if qty_error_param:
            central_quantities[f"{qty_name}_stat_error"] = qty_error_param.GetVal()
    
    # Collect variations by type
    variation_hists = {}
    variation_quantities = {qty: {} for qty in quantity_names}
    
    # Get list of all subdirectories in this mult/pt bin
    parent_dir = input_file.GetDirectory(f"{mult_bin_name}/{pt_bin_name}")
    if not parent_dir:
        print(f"Warning: Parent directory not found for {mult_bin_name}/{pt_bin_name}")
        return
    
    keys = parent_dir.GetListOfKeys()
    
    for key in keys:
        variation_name = key.GetName()
        
        if variation_name == "central":
            continue
        
        var_type = get_variation_type(variation_name)
        if not var_type:
            continue
        
        var_dir = parent_dir.GetDirectory(variation_name)
        if not var_dir:
            continue
        
        # Get variation histogram
        var_hist = var_dir.Get("delta_phi_hist")
        if var_hist:
            var_hist = var_hist.Clone()
            var_hist.SetDirectory(0)
            
            if var_type not in variation_hists:
                variation_hists[var_type] = {}
            variation_hists[var_type][variation_name] = var_hist
        
        # Get variation quantities
        for qty_name in quantity_names:
            qty_param = var_dir.Get(qty_name)
            if qty_param:
                if var_type not in variation_quantities[qty_name]:
                    variation_quantities[qty_name][var_type] = {}
                variation_quantities[qty_name][var_type][variation_name] = qty_param.GetVal()
    
    # Compute systematics for delta phi distribution
    hist_stat, hist_sys, dphi_sys_by_type = compute_deltaphi_systematics(
        central_hist, variation_hists
    )
    
    # Compute systematics for scalar quantities
    scalar_systematics = {}
    scalar_sys_by_type = {}
    
    for qty_name in quantity_names:
        if qty_name in central_quantities:
            sys_unc, sys_by_type = compute_scalar_systematics(
                central_quantities[qty_name],
                variation_quantities[qty_name]
            )
            scalar_systematics[qty_name] = sys_unc
            scalar_sys_by_type[qty_name] = sys_by_type
    
    # Write results to output file
    output_dir_path = f"{mult_bin_name}/{pt_bin_name}"
    output_dir = output_file.mkdir(output_dir_path, output_dir_path, True)
    output_dir.cd()
    
    # Write histograms
    hist_stat.SetName("delta_phi_hist_stat")
    hist_stat.SetTitle(f"#Delta#phi ({mult_bin_name}, {pt_bin_name}) - Statistical Errors")
    hist_stat.Write()
    
    hist_sys.SetName("delta_phi_hist_sys")
    hist_sys.SetTitle(f"#Delta#phi ({mult_bin_name}, {pt_bin_name}) - Systematic Errors")
    hist_sys.Write()
    
    # Write scalar quantities with stat and sys errors
    for qty_name in quantity_names:
        if qty_name in central_quantities:
            # Write value
            param = rt.TParameter(float)(qty_name, central_quantities[qty_name])
            param.Write()
            
            # Write statistical error
            if f"{qty_name}_stat_error" in central_quantities:
                param_stat = rt.TParameter(float)(
                    f"{qty_name}_error_stat",
                    central_quantities[f"{qty_name}_stat_error"]
                )
                param_stat.Write()
            
            # Write systematic error
            if qty_name in scalar_systematics:
                param_sys = rt.TParameter(float)(
                    f"{qty_name}_error_sys",
                    scalar_systematics[qty_name]
                )
                param_sys.Write()
    
    print(f"Processed {mult_bin_name}/{pt_bin_name}")


def main(input_file_path: str = "output/analysis_out.root", 
         output_file_path: str = "output/systematics_out.root"):
    """
    Main function to compute systematic uncertainties.
    
    Args:
        input_file_path: Path to input analysis ROOT file
        output_file_path: Path to output systematics ROOT file
    """
    # Open input and output files
    input_file = rt.TFile(input_file_path, "READ")
    if not input_file or input_file.IsZombie():
        print(f"Error: Could not open input file {input_file_path}")
        return
    
    output_file = rt.TFile(output_file_path, "RECREATE")
    
    # Get list of multiplicity bins
    mult_bins = []
    keys = input_file.GetListOfKeys()
    for key in keys:
        if key.IsFolder():
            mult_bins.append(key.GetName())
    
    # Process each mult/pt bin
    for mult_bin in mult_bins:
        mult_dir = input_file.GetDirectory(mult_bin)
        if not mult_dir:
            continue
        
        # Get pt bins
        pt_keys = mult_dir.GetListOfKeys()
        pt_bins = []
        for key in pt_keys:
            if key.IsFolder():
                pt_bins.append(key.GetName())
        
        for pt_bin in pt_bins:
            process_mult_pt_bin(input_file, output_file, mult_bin, pt_bin)
    
    # Close files
    output_file.Close()
    input_file.Close()
    
    print(f"\nSystematic uncertainties computed and saved to {output_file_path}")


if __name__ == "__main__":
    main()
