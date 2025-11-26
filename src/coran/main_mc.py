from array import array
from collections import defaultdict

import ROOT as rt

from coran import config as cfg

from coran.fitting.mass import FitTypeSignal, FitTypeBg, MassFit
from coran.fitting.correlation import FitTypeJet, FitTypeUE, DeltaPhiFit

from coran.corrections.mixed_event import AcceptanceCorrector, ScaleType
from coran.corrections.sideband import SidebandCorrector

OUTPUT_DICT: defaultdict = defaultdict()
OUTPUT_FILE: rt.TFile | None = None

def write_results_to_root(
    delta_phi_fit: DeltaPhiFit,
    mult_bin_name: str,
    pt_bin_name: str,
    variation_name: str,
):
    """Write analysis results to ROOT file in hierarchical directory structure."""
    global OUTPUT_FILE
    
    if OUTPUT_FILE is None:
        raise RuntimeError("Output file not initialized")
    
    # Create directory structure: mult_bin/pt_bin/variation
    dir_path = f"{mult_bin_name}/{pt_bin_name}/{variation_name}"
    output_dir = OUTPUT_FILE.mkdir(dir_path, dir_path, True)
    output_dir.cd()
    
    # Write the delta phi histogram
    delta_phi_hist = delta_phi_fit.dist.Clone("delta_phi_hist")
    delta_phi_hist.SetTitle(f"#Delta#phi distribution ({mult_bin_name}, {pt_bin_name}, {variation_name})")
    delta_phi_hist.Write()
    
    # Write the fit functions
    fit_total = delta_phi_fit.fit_total.Clone("fit_total")
    fit_total.Write()
    
    fit_jet = delta_phi_fit.fit_jet.Clone("fit_jet")
    fit_jet.Write()
    
    fit_ue = delta_phi_fit.fit_ue.Clone("fit_ue")
    fit_ue.Write()
    
    # Create and write observables as TParameter objects
    # NS width and error
    ns_width = rt.TParameter(float)("ns_width", delta_phi_fit.fit_width_ns.value)
    ns_width.Write()
    ns_width_error = rt.TParameter(float)("ns_width_error", delta_phi_fit.fit_width_ns.stat_error)
    ns_width_error.Write()
    
    # AS width and error
    as_width = rt.TParameter(float)("as_width", delta_phi_fit.fit_width_as.value)
    as_width.Write()
    as_width_error = rt.TParameter(float)("as_width_error", delta_phi_fit.fit_width_as.stat_error)
    as_width_error.Write()
    
    # NS yield and error
    ns_yield = rt.TParameter(float)("ns_yield", delta_phi_fit.ns_yield.value)
    ns_yield.Write()
    ns_yield_error = rt.TParameter(float)("ns_yield_error", delta_phi_fit.ns_yield.stat_error)
    ns_yield_error.Write()
    
    # AS yield and error
    as_yield = rt.TParameter(float)("as_yield", delta_phi_fit.as_yield.value)
    as_yield.Write()
    as_yield_error = rt.TParameter(float)("as_yield_error", delta_phi_fit.as_yield.stat_error)
    as_yield_error.Write()
    
    # UE yield and error
    ue_yield = rt.TParameter(float)("ue_yield", delta_phi_fit.ue_yield.value)
    ue_yield.Write()
    ue_yield_error = rt.TParameter(float)("ue_yield_error", delta_phi_fit.ue_yield.stat_error)
    ue_yield_error.Write()
    
    # Total yield and error
    total_yield = rt.TParameter(float)("total_yield", delta_phi_fit.total_yield.value)
    total_yield.Write()
    total_yield_error = rt.TParameter(float)("total_yield_error", delta_phi_fit.total_yield.stat_error)
    total_yield_error.Write()

def perform_v0_analysis(
    input_correlation_dist: rt.THnSparse,
    input_correlation_dist_mixed: rt.THnSparse,
    num_triggers: float,
    mixed_event_scale_type: ScaleType,
    correlation_fit_type_jet: FitTypeJet,
    correlation_fit_type_ue: FitTypeUE,
    mult_bin_name: str,
    pt_bin_name: str,
    variation_name: str,
):

    signal_min = mass_fit.fit_mean.value + mass_signal_range[0]*mass_fit.fit_width.value
    signal_max = mass_fit.fit_mean.value + mass_signal_range[1]*mass_fit.fit_width.value
    input_correlation_dist.GetAxis(2).SetRangeUser(signal_min, signal_max - cfg.EPSILON)
    input_correlation_dist_mixed.GetAxis(2).SetRangeUser(signal_min, signal_max - cfg.EPSILON)
    signal_3d_same = input_correlation_dist.Projection(0,1,3)
    signal_3d_mixed = input_correlation_dist_mixed.Projection(0,1,3)

    sideband_min = mass_fit.fit_mean.value + mass_sideband_range[0]*mass_fit.fit_width.value
    sideband_max = mass_fit.fit_mean.value + mass_sideband_range[1]*mass_fit.fit_width.value
    input_correlation_dist.GetAxis(2).SetRangeUser(sideband_min, sideband_max - cfg.EPSILON)
    input_correlation_dist_mixed.GetAxis(2).SetRangeUser(sideband_min, sideband_max - cfg.EPSILON)
    sideband_3d = input_correlation_dist.Projection(0,1,3)
    sideband_3d_mixed = input_correlation_dist_mixed.Projection(0,1,3)

    signal_2d_corrector = AcceptanceCorrector(
        signal_3d_same,
        signal_3d_mixed,
        scale_type=ScaleType.ZERO,
        num_bins_z_vertex=cfg.NUM_BINS_Z_VERTEX,
    ) 
    signal_2d_corrector.corrected_2d.Scale(
        1/num_triggers
    )
    signal_2d_corrector.corrected_2d.GetXaxis().SetRangeUser(*cfg.RANGE_DELTA_ETA)
    

    final_dphi = sideband_corrector.corrected_2d.ProjectionY().Clone()

    delta_phi_fit = DeltaPhiFit(
        final_dphi,
        fit_type_jet=correlation_fit_type_jet,
        fit_type_ue=correlation_fit_type_ue,
    )
    
    # Write results to ROOT file
    write_results_to_root(
        delta_phi_fit=delta_phi_fit,
        mult_bin_name=mult_bin_name,
        pt_bin_name=pt_bin_name,
        variation_name=variation_name,
    )
    
    if mult_bin_name not in OUTPUT_DICT:
        OUTPUT_DICT[mult_bin_name] = {}
    if pt_bin_name not in OUTPUT_DICT[mult_bin_name]:
        OUTPUT_DICT[mult_bin_name][pt_bin_name] = {}
    OUTPUT_DICT[mult_bin_name][pt_bin_name][variation_name] = delta_phi_fit
    
def main():
    global OUTPUT_FILE
    
    input_file = rt.TFile(cfg.INPUT_FILE)
    input_list = input_file.Get(cfg.INPUT_LIST)
    
    input_trigger_dist = input_list.FindObject("fTriggerDistEff")
    input_mass_dist = input_list.FindObject("fTriggeredK0Dist")
    
    input_correlation_dist = input_list.FindObject("fDphiHK0Eff")
    input_correlation_dist_mixed = input_list.FindObject("fDphiHK0Mixed")

    OUTPUT_FILE = rt.TFile("output/analysis_out.root", "RECREATE")
    
    for mult_name, (mult_min, mult_max) in cfg.RANGES_MULTIPLICITY.items():
      for pt_name, (pt_min, pt_max) in cfg.RANGES_PT_ASSOCIATED.items():

          input_trigger_dist.GetAxis(cfg.AXIS_HADRON_PT).SetRangeUser(*cfg.RANGE_PT_TRIG)
          input_trigger_dist.GetAxis(cfg.AXIS_HADRON_MULT).SetRangeUser(mult_min, mult_max)
          num_triggers = input_trigger_dist.Projection(cfg.AXIS_HADRON_PT).Integral()

          input_mass_dist.GetAxis(cfg.AXIS_V0_PT).SetRangeUser(pt_min, pt_max)
          input_mass_dist.GetAxis(cfg.AXIS_V0_MULT).SetRangeUser(mult_min, mult_max)
          input_mass_hist = input_mass_dist.Projection(cfg.AXIS_V0_MASS)
          input_mass_hist.SetDirectory(0)

          input_correlation_dist.GetAxis(cfg.AXIS_H_K_PT_TRIGGER).SetRangeUser(*cfg.RANGE_PT_TRIG)
          input_correlation_dist_mixed.GetAxis(cfg.AXIS_H_K_PT_TRIGGER).SetRangeUser(*cfg.RANGE_PT_TRIG)

          input_correlation_dist.GetAxis(cfg.AXIS_H_K_PT_ASSOCIATED).SetRangeUser(pt_min, pt_max)
          input_correlation_dist_mixed.GetAxis(cfg.AXIS_H_K_PT_ASSOCIATED).SetRangeUser(pt_min, pt_max)

          input_correlation_dist.GetAxis(cfg.AXIS_H_K_MULT).SetRangeUser(mult_min, mult_max)
          input_correlation_dist_mixed.GetAxis(cfg.AXIS_H_K_MULT).SetRangeUser(mult_min, mult_max)

          axes = array(
              "i",
              [
                2, 3, 4, 5
              ],
          )
          input_correlation_dist_proj = input_correlation_dist.Projection(len(axes), axes)
          input_correlation_dist_mixed_proj = input_correlation_dist_mixed.Projection(len(axes), axes)

          perform_v0_analysis(
              input_mass_hist = input_mass_hist,
              input_correlation_dist = input_correlation_dist_proj,
              input_correlation_dist_mixed = input_correlation_dist_mixed_proj,
              particle_type="k0",
              num_triggers=num_triggers,
              mass_fit_type_signal=FitTypeSignal.CRYSTAL,
              mass_fit_type_bg=FitTypeBg.POL1,
              mass_signal_range=cfg.RANGE_SIGNAL,
              mass_sideband_range=cfg.RANGE_SIDEBAND,
              mixed_event_scale_type=ScaleType.ZERO,
              correlation_fit_type_jet=FitTypeJet.VON,
              correlation_fit_type_ue=FitTypeUE.AVG_4,
              mult_bin_name=mult_name,
              pt_bin_name=pt_name,
              variation_name="central"
          )
          
          for variation_name, mass_signal_range in cfg.VARIATIONS_SIGNAL.items():
              perform_v0_analysis(
                  input_mass_hist = input_mass_hist,
                  input_correlation_dist = input_correlation_dist_proj,
                  input_correlation_dist_mixed = input_correlation_dist_mixed_proj,
                  particle_type="k0",
                  num_triggers=num_triggers,
                  mass_fit_type_signal=FitTypeSignal.CRYSTAL,
                  mass_fit_type_bg=FitTypeBg.POL1,
                  mass_signal_range=mass_signal_range,
                  mass_sideband_range=cfg.RANGE_SIDEBAND,
                  mixed_event_scale_type=ScaleType.ZERO,
                  correlation_fit_type_jet=FitTypeJet.VON,
                  correlation_fit_type_ue=FitTypeUE.AVG_4,
                  mult_bin_name=mult_name,
                  pt_bin_name=pt_name,
                  variation_name="variation_signal_" + variation_name
              )

          for variation_name, mass_sideband_range in cfg.VARIATIONS_SIDEBAND.items():
                 perform_v0_analysis(
                    input_mass_hist = input_mass_hist,
                    input_correlation_dist = input_correlation_dist_proj,
                    input_correlation_dist_mixed = input_correlation_dist_mixed_proj,
                    particle_type="k0",
                    num_triggers=num_triggers,
                    mass_fit_type_signal=FitTypeSignal.CRYSTAL,
                    mass_fit_type_bg=FitTypeBg.POL1,
                    mass_signal_range=cfg.RANGE_SIGNAL,
                    mass_sideband_range=mass_sideband_range,
                    mixed_event_scale_type=ScaleType.ZERO,
                    correlation_fit_type_jet=FitTypeJet.VON,
                    correlation_fit_type_ue=FitTypeUE.AVG_4,
                    mult_bin_name=mult_name,
                    pt_bin_name=pt_name,
                    variation_name="variation_sideband_" + variation_name
                )
  
    OUTPUT_FILE.Close()
    
    print(f"Analysis complete. Results written to output/analysis_out.root")
    
if __name__ == "__main__":
    main()
