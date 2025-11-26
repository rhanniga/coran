import ROOT as rt

rt.gROOT.SetBatch(True)

from coran.models import StartingPar
from coran.fitting.correlation import FitTypeJet, FitTypeUE, DeltaPhiFit

# --- User input: choose plot mode ---
print("Choose plotting mode:")
print("1. Plot NS and AS separately")
print("2. Plot combined Jet-Like only")
plot_mode = input("Enter your choice (1 or 2): ").strip()

plot_separate = (plot_mode == "1")

# --- 1. Load histograms ---
f = rt.TFile("phsd_dphi_new.root")

# Define multiplicity bins
mult_bins = ["0_20", "20_50", "50_80"]
mult_labels = ["0-20%", "20-50%", "50-80%"]

# Load histograms for each multiplicity bin
h_h_dphi_dict = {}
h_k0_dphi_dict = {}

for mult_bin in mult_bins:
    h_h_name = f"h_h_dphi_phsd_cent_{mult_bin}_trigger_4_8_assoc_1_4_delta_eta_12"
    h_k0_name = f"h_k0_dphi_phsd_cent_{mult_bin}_trigger_4_8_assoc_1_4_delta_eta_12"
    
    h_h_dphi_dict[mult_bin] = f.Get(h_h_name)
    h_k0_dphi_dict[mult_bin] = f.Get(h_k0_name)
    
    h_h_dphi_dict[mult_bin].SetDirectory(0)
    h_k0_dphi_dict[mult_bin].SetDirectory(0)

f.Close()

# --- 2. Define starting parameters and perform fits ---
# For Von Mises + V2 UE fit, we have 7 parameters:
# 0: NS amplitude (yield)
# 1: NS kappa (width)
# 2: AS amplitude (yield)
# 3: AS kappa (width)
# 4: UE pedestal
# 5: v2_trig
# 6: v2_assoc

fit_type_jet = FitTypeJet.VON
fit_type_ue = FitTypeUE.V2

# Perform fits for each multiplicity bin
fits_h_h = {}
fits_h_k0 = {}

for mult_bin in mult_bins:
    print(f"\nFitting {mult_bin} multiplicity bin...")
    print(f"Fitting h-h correlation...")
    fits_h_h[mult_bin] = DeltaPhiFit(h_h_dphi_dict[mult_bin], fit_type_jet, fit_type_ue, (0.05, 0.05))
    print(f"Fitting h-K0s correlation...")
    fits_h_k0[mult_bin] = DeltaPhiFit(h_k0_dphi_dict[mult_bin], fit_type_jet, fit_type_ue, (0.05, 0.05))


# --- 3. Extract yields and calculate per-trigger yields ---
def get_yields(fit: DeltaPhiFit, hist: rt.TH1D):
    """Extracts yields for different regions."""
    n_trig = 1 # Assuming 1 trigger for PHSD, adjust if not the case
    
    # Jet yields are the amplitudes from the fit
    ns_yield = fit.fit_amplitude_ns
    as_yield = fit.fit_amplitude_as
    
    jet_like_yield = ns_yield.value + as_yield.value
    jet_like_error = rt.TMath.Sqrt(ns_yield.stat_error**2 + as_yield.stat_error**2)

    # UE yield is the integral of the UE function
    # For V2, UE has parameters: pedestal, v2_trig, v2_assoc
    # But fit.fit_ue is just the UE component function, so parameter 0 should be the pedestal
    # However, we should get it from the total fit function at parameter index 4
    ue_fit_func = fit.fit_total  # Use total fit function
    pedestal = ue_fit_func.GetParameter(4)  # UE pedestal is parameter 4
    pedestal_error = ue_fit_func.GetParError(4)
    
    # Scale up UE error for visibility (artificially)
    ue_error_scale_factor = 5.0  # Increase this to make UE band more visible
    pedestal_error = pedestal_error * ue_error_scale_factor
    
    # The range of the dphi distribution is 2*pi
    phi_range = hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()
    
    ue_yield = pedestal * phi_range
    ue_yield_error = pedestal_error * phi_range
    
    total_yield = ns_yield.value + as_yield.value + ue_yield
    total_yield_error = rt.TMath.Sqrt(ns_yield.stat_error**2 + as_yield.stat_error**2 + ue_yield_error**2)

    yields = {
        "UE": {"yield": ue_yield / n_trig, "error": ue_yield_error / n_trig},
        "Total": {"yield": total_yield / n_trig, "error": total_yield_error / n_trig}
    }
    
    if plot_separate:
        yields["NS"] = {"yield": ns_yield.value / n_trig, "error": ns_yield.stat_error / n_trig}
        yields["AS"] = {"yield": as_yield.value / n_trig, "error": as_yield.stat_error / n_trig}
    else:
        yields["Jet-Like"] = {"yield": jet_like_yield / n_trig, "error": jet_like_error / n_trig}
    
    return yields

# Calculate yields and ratios for each multiplicity bin
yields_h_h = {}
yields_h_k0 = {}
ratios_by_mult = {}

for mult_bin in mult_bins:
    yields_h_h[mult_bin] = get_yields(fits_h_h[mult_bin], h_h_dphi_dict[mult_bin])
    yields_h_k0[mult_bin] = get_yields(fits_h_k0[mult_bin], h_k0_dphi_dict[mult_bin])
    
    # Calculate ratios for this multiplicity bin
    if plot_separate:
        regions = ["NS", "AS", "UE", "Total"]
    else:
        regions = ["Jet-Like", "UE", "Total"]
    
    ratios_by_mult[mult_bin] = {}
    
    for region in regions:
        y_k0 = yields_h_k0[mult_bin][region]["yield"]
        e_k0 = yields_h_k0[mult_bin][region]["error"]
        y_hh = yields_h_h[mult_bin][region]["yield"]
        e_hh = yields_h_h[mult_bin][region]["error"]
        
        # Debug: print individual errors
        if region == "UE":
            print(f"{mult_bin} UE - y_k0: {y_k0:.6f}, e_k0: {e_k0:.6f}, y_hh: {y_hh:.6f}, e_hh: {e_hh:.6f}")
        
        ratio = y_k0 / y_hh if y_hh != 0 else 0
        # Calculate error with proper handling of zero errors
        if y_k0 != 0 and y_hh != 0 and (e_k0 != 0 or e_hh != 0):
            error = ratio * rt.TMath.Sqrt((e_k0/y_k0)**2 + (e_hh/y_hh)**2)
        else:
            error = 0
        ratios_by_mult[mult_bin][region] = {"ratio": ratio, "error": error}
        print(f"{mult_bin} {region} Ratio (h-K0s/h-h): {ratio:.3f} +/- {error:.3f}")

# --- 5. Plotting ---
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptFit(0)

# Define plot correlation function
def plot_correlation(canvas, hist, fit, title, color, label, mult_label):
    canvas.cd()
    hist.SetTitle("")
    hist.GetXaxis().SetTitle("#Delta#varphi")
    hist.GetXaxis().SetTitleOffset(0.9)
    hist.GetYaxis().SetTitle("1/#it{N}_{trig} d#it{N}/d#Delta#phi")
    hist.GetYaxis().SetTitleOffset(1.4)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1.5)
    hist.SetMarkerColor(rt.kBlack)
    hist.SetLineColor(rt.kBlack)
    hist.GetYaxis().SetRangeUser(hist.GetMinimum()*0.9, hist.GetMaximum() * 1.1)
    hist.Draw("PE")
    
    pt_label = rt.TLatex()
    pt_label.SetNDC()
    pt_label.SetTextSize(0.045)
    pt_label.DrawLatex(0.54, 0.6, "#bf{4.0 < #it{p}_{T,trig} < 8.0 GeV/#it{c}}")
    pt_label.DrawLatex(0.54, 0.52, "#bf{1.0 < #it{p}_{T,assoc} < 4.0 GeV/#it{c}}")

    phsd_label = rt.TLatex()
    phsd_label.SetNDC()
    phsd_label.SetTextSize(0.06)
    phsd_label.DrawLatex(0.18, 0.82, f"PHSD {label}")
    
    mult_text = rt.TLatex()
    mult_text.SetNDC()
    mult_text.SetTextSize(0.05)
    mult_text.DrawLatex(0.18, 0.74, f"Mult. Percentile: {mult_label}")
    
    fit.fit_total.SetLineColor(color)
    fit.fit_total.SetLineStyle(rt.kSolid)
    fit.fit_total.Draw("SAME")
    
    fit.fit_ue.SetLineColor(color)
    fit.fit_ue.SetLineStyle(rt.kDashed)
    fit.fit_ue.Draw("SAME")
    
    legend = rt.TLegend(0.6, 0.7, 0.9, 0.9)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.AddEntry(hist, label, "lep")
    legend.AddEntry(fit.fit_total, "Total fit", "l")
    legend.AddEntry(fit.fit_ue, "UE fit", "l")
    legend.Draw()
    canvas.Update()
    canvas.SaveAs(f"{title.replace(' ', '_').replace('%', '')}.pdf")
    return legend, pt_label, phsd_label, mult_text # Keep objects in scope

# Plot correlations for each multiplicity bin
for i, mult_bin in enumerate(mult_bins):
    mult_label = mult_labels[i]
    
    # h-h correlation
    h_h_canvas = rt.TCanvas(f"h_h_canvas_{mult_bin}", "", 1600, 1200)
    legend_hh, pt_hh, phsd_hh, mult_hh = plot_correlation(
        h_h_canvas, h_h_dphi_dict[mult_bin], fits_h_h[mult_bin], 
        f"h-h Correlation {mult_label}", rt.kRed + 1, "h-h", mult_label
    )
    
    # h-K0s correlation
    h_k0_canvas = rt.TCanvas(f"h_k0_canvas_{mult_bin}", "", 1600, 1200)
    legend_k0, pt_k0, phsd_k0, mult_k0 = plot_correlation(
        h_k0_canvas, h_k0_dphi_dict[mult_bin], fits_h_k0[mult_bin], 
        f"h-K Correlation {mult_label}", rt.kBlue + 1, "h-K", mult_label
    )

# Plot ratios vs multiplicity
def plot_ratios_vs_mult(canvas, ratios_by_mult, mult_bins, mult_labels, regions):
    """Plots the yield ratios vs multiplicity on a canvas."""
    canvas.cd()
    
    # Use TGraphErrors for plotting
    graphs = {}
    if plot_separate:
        colors = {"NS": rt.kRed + 1, "AS": rt.kBlue + 1, "UE": rt.kGreen + 2}
        markers = {"NS": 20, "AS": 21, "UE": 22}
        legend_labels = {"NS": "NS", "AS": "AS", "UE": "UE"}
    else:
        colors = {"Jet-Like": rt.kViolet + 1, "UE": rt.kGreen + 2}
        markers = {"Jet-Like": 20, "UE": 22}
        legend_labels = {"Jet-Like": "Jet", "UE": "UE"}

    # Create a dummy histogram for axis (reversed order)
    n_mult_bins = len(mult_bins)
    h_axis = rt.TH1F("h_axis", ";Mult. Percentile;(h-K)/(h-h) yield ratio", n_mult_bins, 0, n_mult_bins)
    
    # Set bin labels in reverse order
    for i, mult_label in enumerate(reversed(mult_labels)):
        h_axis.GetXaxis().SetBinLabel(i + 1, mult_label)
    
    # Find y-axis range
    all_ratios = []
    for mult_bin in mult_bins:
        for region in regions:
            if region != "Total":
                all_ratios.append(ratios_by_mult[mult_bin][region]["ratio"])
    
    y_min = min(all_ratios) * 0.8 if all_ratios else 0
    y_max = max(all_ratios) * 1.2 if all_ratios else 1
    h_axis.SetMinimum(y_min)
    h_axis.SetMaximum(y_max)
    h_axis.GetYaxis().SetTitleOffset(1.4)
    h_axis.Draw("axis")

    # Create graphs for each region
    for region in regions:
        if region == "Total":
            continue
            
        graph = rt.TGraphErrors(n_mult_bins)
        # Reverse the order when plotting
        for i, mult_bin in enumerate(reversed(mult_bins)):
            ratio = ratios_by_mult[mult_bin][region]["ratio"]
            error = ratios_by_mult[mult_bin][region]["error"]
            print(region, error)
            graph.SetPoint(i, i + 0.5, ratio)
            graph.SetPointError(i, 0.4, error)  # Set x-error to create horizontal band
        
        graph.SetMarkerColor(colors.get(region, rt.kBlack))
        graph.SetLineColor(colors.get(region, rt.kBlack))
        graph.SetMarkerStyle(markers.get(region, 20))
        graph.SetMarkerSize(2.0)
        graph.SetLineWidth(2)
        graph.SetFillColorAlpha(colors.get(region, rt.kBlack), 0.3)  # Increased opacity
        graph.Draw("P E3 SAME")  # E4 draws a smoothed filled band
        graphs[region] = graph

    # Add pt labels
    pt_label = rt.TLatex()
    pt_label.SetNDC()
    pt_label.SetTextSize(0.060)
    pt_label.DrawLatex(0.18, 0.83, "PHSD (h-K)/(h-h)")
    pt_label.SetTextSize(0.045)
    pt_label.DrawLatex(0.18, 0.75, "#bf{4.0 < #it{p}_{T,trig} < 8.0 GeV/#it{c}}")
    pt_label.DrawLatex(0.18, 0.68, "#bf{1.0 < #it{p}_{T,assoc} < 4.0 GeV/#it{c}}")

    legend = rt.TLegend(0.70, 0.68, 0.9, 0.85)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    for region, graph in reversed(graphs.items()):
        legend.AddEntry(graph, legend_labels[region], "pf")
    legend.Draw()

    canvas.Update()
    canvas.SaveAs("yield_ratios_vs_multiplicity.pdf")
    return graphs, legend, h_axis, pt_label  # Keep objects in scope

canvas_ratios_mult = rt.TCanvas("canvas_ratios_mult", "Yield Ratios vs Multiplicity", 1000, 800)
if plot_separate:
    regions = ["NS", "AS", "UE", "Total"]
else:
    regions = ["Jet-Like", "UE", "Total"]
graphs_mult, legend_mult, axis_mult, pt_mult = plot_ratios_vs_mult(canvas_ratios_mult, ratios_by_mult, mult_bins, mult_labels, regions)

print("\nPlots saved successfully!")
print("- Delta-phi correlations for each multiplicity bin")
print("- yield_ratios_vs_multiplicity.pdf")

