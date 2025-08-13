from uuid import uuid4, UUID
from enum import Enum

import ROOT as rt

from coran.fitting.common import StartingPar, Observable

# TODO: think about where to put this
BINS_AVG_4 = [1, 8, 9, 16]
BINS_AVG_6 = [1, 2, 7, 8, 9, 16]


class FitTypeUE(Enum):
    V2 = "v_{2}"
    AVG_4 = "4 bin avg"
    AVG_6 = "6 bin avg"


class FitTypeJet(Enum):
    GAUS = "Gaussian"
    VON = "Von Mises"


def get_width_from_kappa(kappa):
    return rt.TMath.Sqrt(
        -2 * rt.TMath.Log(rt.TMath.BesselI1(kappa) / rt.TMath.BesselI0(kappa))
    )


def get_width_error_from_kappa(kappa, kappa_error):
    deriv = (
        rt.TMath.BesselI0(kappa) ** 2
        + rt.TMath.BesselI0(kappa) * rt.TMath.BesselI(2, kappa)
        - 2 * rt.TMath.BesselI1(kappa) ** 2
    ) / (
        2
        * rt.TMath.Sqrt(2)
        * rt.TMath.BesselI0(kappa)
        * rt.TMath.BesselI1(kappa)
        * rt.TMath.Sqrt(
            -rt.TMath.Log(rt.TMath.BesselI1(kappa) / rt.TMath.BesselI0(kappa))
        )
    )
    return deriv * kappa_error


def get_ue_avg(
    delta_phi_dist: rt.TH1D,
    num_bins: int = 6,
) -> tuple[float, float]:
    assert num_bins in [4, 6], "only 4 or 6 bins are supported for average bg"
    if num_bins == 4:
        ue_avg = sum(
            delta_phi_dist.GetBinContent(bin_index) for bin_index in BINS_AVG_4
        )
        ue_avg /= 4

        ue_avg_error = sum(
            delta_phi_dist.GetBinError(bin_index) ** 2 for bin_index in BINS_AVG_4
        )
        ue_avg_error = (1 / 4) * rt.TMath.Sqrt(ue_avg_error)

    else:
        ue_avg = sum(
            delta_phi_dist.GetBinContent(bin_index) for bin_index in BINS_AVG_6
        )
        ue_avg /= 6

        ue_avg_error = sum(
            delta_phi_dist.GetBinError(bin_index) ** 2 for bin_index in BINS_AVG_6
        )
        ue_avg_error = (1 / 6) * rt.TMath.Sqrt(ue_avg_error)

    return ue_avg, ue_avg_error


class FitDescriptor:
    def __init__(self, fit_type_jet: FitTypeJet, fit_type_ue: FitTypeUE):

        self._fit_type_jet = fit_type_jet
        self._fit_type_ue = fit_type_ue

        self._num_pars = 0
        self._par_index_start_jet = 0
        self._par_index_start_ue = 0

        self._par_index_amplitude_ns = 0
        self._par_index_width_ns = 0

        self._par_index_amplitude_as = 0
        self._par_index_width_as = 0

        self._generate_fit_strings()
        self._generate_labels()

    def _generate_fit_strings(self):
        fit_par_index = 0
        self._par_index_start_jet = fit_par_index
        if self._fit_type_jet == FitTypeJet.GAUS:
            self._fit_string_jet = (
                f"gaus({fit_par_index}) + gaus({fit_par_index + 3}) + "
                f"gaus({fit_par_index + 6}) + gaus({fit_par_index + 9})"
            )
            fit_par_index += 12

            # not sure if there's anything to do other than to hard-code these...
            self._par_index_amplitude_ns = 0
            self._par_index_width_ns = 1

            self._par_index_amplitude_as = 6
            self._par_index_width_as = 7

        elif self._fit_type_jet == FitTypeJet.VON:
            self._fit_string_jet = (
                f"[{fit_par_index}]/(2*TMath::Pi()*TMath::BesselI0([{fit_par_index + 1}]))*"
                f"TMath::Exp([{fit_par_index + 1}]*TMath::Cos(x)) + "
                f"[{fit_par_index + 2}]/(2*TMath::Pi()*TMath::BesselI0([{fit_par_index + 3}]))*"
                f"TMath::Exp([{fit_par_index + 3}]*TMath::Cos(x - TMath::Pi()))"
            )
            fit_par_index += 4

            # not sure if there's anything to do other than to hard-code these...
            self._par_index_amplitude_ns = 0
            self._par_index_width_ns = 1

            self._par_index_amplitude_as = 2
            self._par_index_width_as = 3

        self._fit_string_total = self._fit_string_jet

        self._par_index_start_ue = fit_par_index

        if self._fit_type_ue == FitTypeUE.V2:
            self._fit_string_total += f"+ [{fit_par_index}]*(1 + 2*([{fit_par_index + 1}]*[{fit_par_index + 2}]*cos(2*x)))"
            self._fit_string_ue = f"[0]*(1 + 2*([1]*[2]*cos(2*x)))"
            fit_par_index += 3
        elif self._fit_type_ue in [FitTypeUE.AVG_4, FitTypeUE.AVG_6]:
            self._fit_string_total += f"+ [{fit_par_index}]"
            self._fit_string_ue = "[0]"
            fit_par_index += 1

        self._num_pars = fit_par_index

    def _generate_labels(self):
        self._fit_label_total = (
            f"{self._fit_type_jet.value} + {self._fit_type_ue.value}"
        )
        self._fit_label_jet = self._fit_type_jet.value
        self._fit_label_ue = self._fit_type_ue.value

    @property
    def fit_string_total(self):
        return self._fit_string_total

    @property
    def fit_string_jet(self):
        return self._fit_string_jet

    @property
    def fit_string_ue(self):
        return self._fit_string_ue

    @property
    def fit_label_total(self):
        return self._fit_label_total

    @property
    def fit_label_jet(self):
        return self._fit_label_jet

    @property
    def fit_label_ue(self):
        return self._fit_label_ue

    @property
    def num_pars(self) -> int:
        return self._num_pars

    @property
    def par_index_start_jet(self) -> int:
        return self._par_index_start_jet

    @property
    def par_index_start_ue(self) -> int:
        return self._par_index_start_ue

    @property
    def par_index_amplitude_ns(self) -> int:
        return self._par_index_amplitude_ns

    @property
    def par_index_width_ns(self) -> int:
        return self._par_index_width_ns

    @property
    def par_index_amplitude_as(self) -> int:
        return self._par_index_amplitude_as

    @property
    def par_index_width_as(self) -> int:
        return self._par_index_width_as


class DeltaPhiFit:
    def __init__(
        self,
        dist: rt.TH1D,
        fit_type_jet: FitTypeJet,
        fit_type_ue: FitTypeUE,
        starting_pars: list[StartingPar],
    ) -> None:

        self._id = uuid4()
        self._dist = dist.Clone(f"{self._id}_dist")
        self._fit_type_jet = fit_type_jet
        self._fit_type_ue = fit_type_ue
        self._starting_pars = starting_pars

        self._fit_descriptor = FitDescriptor(
            fit_type_jet=self._fit_type_jet, fit_type_ue=self._fit_type_ue
        )

        assert self._fit_descriptor.num_pars == len(self._starting_pars), (
            f"Number of starting parameters ({len(self._starting_pars)}) does not match "
            f"the number of parameters in the fit descriptor ({self._fit_descriptor.num_pars})."
        )

        self._ue_avg, self._ue_error = get_ue_avg(
            self._dist, num_bins=4 if fit_type_ue == FitTypeUE.AVG_4 else 6
        )

        self._fit_total = None
        self._fit_ue = None
        self._fit_jet = None

        self._configure_fits()

        self._fit_dist()

    def _configure_fits(self):
        self._configure_fit_total()
        self._configure_fit_ue()
        self._configure_fit_jet()

    def _configure_fit_total(self):
        fit_total = rt.TF1(
            f"{self._id}_fit_total",
            self._fit_descriptor.fit_string_total,
            self._dist.GetXaxis().GetXmin(),
            self._dist.GetXaxis().GetXmax(),
        )
        for par_index, starting_par in enumerate(self._starting_pars):
            if starting_par.limits is not None:
                fit_total.SetParLimits(par_index, *starting_par.limits)
            if starting_par.fixed:
                fit_total.FixParameter(par_index, starting_par.value)
            else:
                fit_total.SetParameter(par_index, starting_par.value)

        # for now, we always set the pedestal to the average before fitting
        fit_total.SetParameter(self._fit_descriptor.par_index_start_ue, self._ue_avg)
        fit_total.SetParLimits(
            self._fit_descriptor.par_index_start_ue,
            self._ue_avg - 2 * self._ue_error,
            self._ue_avg + 2 * self._ue_error,
        )

        self._fit_total = fit_total
        self._fit_total.SetNpx(1000)

    def _configure_fit_ue(self):
        self._fit_ue = rt.TF1(
            f"{self._id}_fit_ue",
            self._fit_descriptor.fit_string_ue,
            self._dist.GetXaxis().GetXmin(),
            self._dist.GetXaxis().GetXmax(),
        )
        self._fit_ue.SetNpx(1000)

    def _configure_fit_jet(self):
        self._fit_jet = rt.TF1(
            f"{self._id}_fit_jet",
            self._fit_descriptor.fit_string_jet,
            self._dist.GetXaxis().GetXmin(),
            self._dist.GetXaxis().GetXmax(),
        )
        self._fit_jet.SetNpx(1000)

    def _fit_dist(self):
        self._fit_result = self._dist.Fit(
            self._fit_total,
            "RBS",
            "",
            self._dist.GetXaxis().GetXmin(),
            self._dist.GetXaxis().GetXmax(),
        )

        for par_index in range(
            self._fit_descriptor.par_index_start_jet,
            self._fit_descriptor.par_index_start_ue,
        ):
            self._fit_jet.SetParameter(
                par_index - self._fit_descriptor.par_index_start_jet,
                self._fit_total.GetParameter(par_index),
            )

        for par_index in range(
            self._fit_descriptor.par_index_start_ue, self._fit_descriptor.num_pars
        ):
            self._fit_ue.SetParameter(
                par_index - self._fit_descriptor.par_index_start_ue,
                self._fit_total.GetParameter(par_index),
            )

    @property
    def fit_width_ns(self) -> Observable:
        width_ns = self._fit_total.GetParameter(self._fit_descriptor.par_index_width_ns)

        width_ns_error = self._fit_total.GetParError(
            self._fit_descriptor.par_index_width_ns
        )

        return Observable(
            value=(
                width_ns
                if self._fit_type_jet != FitTypeJet.VON
                else get_width_from_kappa(width_ns)
            ),
            stat_error=(
                width_ns_error
                if self._fit_type_jet != FitTypeJet.VON
                else get_width_error_from_kappa(width_ns, width_ns_error)
            ),
        )

    @property
    def fit_amplitude_ns(self) -> Observable:
        ns_amplitude = self._fit_total.GetParameter(
            self._fit_descriptor.par_index_amplitude_ns
        )
        ns_amplitude_error = self._fit_total.GetParError(
            self._fit_descriptor.par_index_amplitude_ns
        )
        return Observable(value=ns_amplitude, stat_error=ns_amplitude_error)

    @property
    def fit_width_as(self) -> Observable:
        width_as = self._fit_total.GetParameter(self._fit_descriptor.par_index_width_as)
        width_as_error = self._fit_total.GetParError(
            self._fit_descriptor.par_index_width_as
        )

        return Observable(
            value=(
                width_as
                if self._fit_type_jet != FitTypeJet.VON
                else get_width_from_kappa(width_as)
            ),
            stat_error=(
                width_as_error
                if self._fit_type_jet != FitTypeJet.VON
                else get_width_error_from_kappa(width_as, width_as_error)
            ),
        )

    @property
    def fit_amplitude_as(self) -> Observable:
        as_amplitude = self._fit_total.GetParameter(
            self._fit_descriptor.par_index_amplitude_as
        )
        as_amplitude_error = self._fit_total.GetParError(
            self._fit_descriptor.par_index_amplitude_as
        )
        return Observable(value=as_amplitude, stat_error=as_amplitude_error)

    @property
    def fit_total(self) -> rt.TF1:
        return self._fit_total

    @property
    def fit_ue(self) -> rt.TF1:
        return self._fit_ue

    @property
    def fit_jet(self) -> rt.TF1:
        return self._fit_jet

    @property
    def fit_label_total(self) -> str:
        return self._fit_descriptor.fit_label_total

    @property
    def fit_label_ue(self) -> str:
        return self._fit_descriptor.fit_label_ue

    @property
    def fit_label_jet(self) -> str:
        return self._fit_descriptor.fit_label_jet
