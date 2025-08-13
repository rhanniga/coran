import ROOT as rt

rt.gROOT.SetBatch(True)

import coran.config as cfg

from coran.corrections.mixed_event import AcceptanceCorrector
from coran.corrections.sideband import SidebandCorrector

from coran.ranges import RangeMassK0

from coran.fitting.common import StartingPar, Observable
from coran.fitting.correlation import DeltaPhiFit
from coran.fitting.mass import MassFit, FitTypeSignal, FitTypeBg


input_file = rt.TFile("input/" + cfg.INPUT_FILENAME, "READ")
input_list = input_file.Get(cfg.INPUT_LIST) if cfg.INPUT_LIST else input_file

trigger_dist = input_list.FindObject("fTriggerDistEff")

k0_dist = input_list.FindObject("fTriggeredK0Dist")
h_k0_dist = input_list.FindObject("fDphiHK0Eff")
h_k0_mixed_dist = input_list.FindObject("fDphiHK0Mixed")


h_k0_dist.GetAxis(cfg.AXIS_H_K_PT_TRIGGER).SetRangeUser(4.0, 7.999)
h_k0_mixed_dist.GetAxis(cfg.AXIS_H_K_PT_TRIGGER).SetRangeUser(4.0, 7.999)

h_k0_dist_6d = h_k0_dist.Projection(1, 2, 3, 4, 5, 6)

# hard code the axes to test some things, this is definitely nowhere near finalized
for range_mult in cfg.RANGES_MULTIPLICITY:
    trigger_dist.GetAxis(cfg.AXIS_HADRON_MULT).SetRangeUser(*range_mult)
    k0_dist.GetAxis(cfg.AXIS_V0_MULT).SetRangeUser(*range_mult)

    h_k0_dist.GetAxis(cfg.AXIS_H_K_MULT).SetRangeUser(*range_mult)
    h_k0_mixed_dist.GetAxis(cfg.AXIS_H_K_MULT).SetRangeUser(*range_mult)

    h_k0_dist_6d.GetAxis(5).SetRangeUser(*range_mult)

    h_k0_dist_5d = h_k0_dist_6d.Projection(0, 1, 2, 3, 4)

    for range_pt_assoc in cfg.RANGES_PT_ASSOCIATED:
        h_k0_dist.GetAxis(cfg.AXIS_H_K_PT_ASSOCIATED).SetRangeUser(*range_pt_assoc)
        h_k0_mixed_dist.GetAxis(cfg.AXIS_H_K_PT_ASSOCIATED).SetRangeUser(
            *range_pt_assoc
        )
        k0_dist.GetAxis(cfg.AXIS_V0_PT).SetRangeUser(*range_pt_assoc)

        h_k0_dist_5d.GetAxis(0).SetRangeUser(*range_pt_assoc)
        h_k0_dist_4d = h_k0_dist_5d.Projection(1, 2, 3, 4)

        # h_k0_mass = h_k0_dist.Projection(cfg.AXIS_H_K_MASS, "O")
        h_k0_mass = h_k0_dist_4d.Projection(2)
        # h_k0_mass = k0_dist.Projection(cfg.AXIS_V0_MASS)

        h_k0_mass_fit = MassFit(
            hist=h_k0_mass,
            fit_type_signal=FitTypeSignal.CRYSTAL,
            fit_type_bg=FitTypeBg.POL1,
            starting_pars=[
                StartingPar(
                    h_k0_mass.GetMaximum(),
                    (1.0 * h_k0_mass.GetMaximum(), 1.5 * h_k0_mass.GetMaximum()),
                ),  # yield
                StartingPar(0.4976, (0.49, 0.50)),  # mean
                StartingPar(0.006),  # sigma
                StartingPar(1.2),  # alpha_l
                StartingPar(0.9),  # alpha_h
                StartingPar(6.0),  # n_l
                StartingPar(10.0),  # n_h
                StartingPar(0.0),
                StartingPar(0.0),
            ],
            range_signal=RangeMassK0(-5.0, 5.0),
            range_sideband=RangeMassK0(-16.5, -6.5),
            fit_range=(0.44, 0.56 - cfg.EPSILON),
        )

        canvas = rt.TCanvas(f"{range_pt_assoc}_{range_mult}_canvas", "", 800, 600)

        h_k0_mass.Draw()
        h_k0_mass_fit.fit_total.Draw("same")

        canvas.SaveAs(f"output/{range_pt_assoc}_{range_mult}_mass_fit.png")
