import ROOT as rt

rt.gROOT.SetBatch(True)

from coran.fitting.mass import FitTypeSignal, FitTypeBg, MassFit


input_file = rt.TFile("test_mass.root")
input_dist = input_file.Get("test_mass_dist")

test_fit = MassFit(
    input_dist,
    particle_type="k0",
    range_signal=(0.48, 0.515),
    range_sideband=(0.44, 0.465),
    fit_type_signal=FitTypeSignal.CRYSTAL,
    fit_type_bg=FitTypeBg.POL1,
    debug=True,
)


