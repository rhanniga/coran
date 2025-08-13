from dataclasses import dataclass
from math import exp
from enum import Enum
from uuid import uuid4, UUID

import ROOT as rt

from coran.ranges import RangeMassK0
from coran.fitting.common import Observable
from coran.fitting.common import StartingPar

FUN_PARS_CRYSTAL = 7
FUN_PARS_GAUS = 3


class FitTypeSignal(Enum):
    GAUS = "Gaussian"
    CRYSTAL = "Double Crystal Ball"


class FitTypeBg(Enum):
    POL1 = "Linear"
    POL2 = "Quadratic"


@rt.Numba.Declare(["double"] * FUN_PARS_CRYSTAL, "double")
def double_crystal(
    x: float,
    mean: float,
    sigma: float,
    alpha_l: float,
    alpha_h: float,
    n_l: float,
    n_h: float,
) -> float:

    if sigma == 0 or n_l == 0 or n_h == 0 or alpha_l == 0 or alpha_h == 0:
        return 0.0

    arg = (x - mean) / sigma
    f1_l, f1_h = alpha_l / n_l, alpha_h / n_h
    f2_l = (n_l / alpha_l) - alpha_l - arg
    f2_h = (n_h / alpha_h) - alpha_h + arg

    if -alpha_l <= arg <= alpha_h:
        return exp(-0.5 * arg**2)
    if arg < -alpha_l:
        return exp(-0.5 * alpha_l**2) * pow(f1_l * f2_l, -n_l)
    return exp(-0.5 * alpha_h**2) * pow(f1_h * f2_h, -n_h)


@rt.Numba.Declare(["double"] * FUN_PARS_GAUS, "double")
def gaus(x: float, mean: float, sigma: float) -> float:

    if sigma == 0:
        return 1.0e30

    arg = (x - mean) / sigma

    if abs(arg) > 39.0:
        return 0.0
    return exp(-0.5 * arg * arg)


def get_straight_line(point1, point2):
    # Ensure the points are valid and not identical
    if point1[0] == point2[0]:
        raise ValueError("The x-coordinates of the two points cannot be the same.")

    # Calculate the slope (m) and intercept (b)
    m = (point2[1] - point1[1]) / (point2[0] - point1[0])
    b = point1[1] - m * point1[0]

    return b, m


class FitDescriptor:
    def __init__(self, fit_type_signal: FitTypeSignal, fit_type_bg: FitTypeBg):
        self._fit_type_signal = fit_type_signal
        self._fit_type_bg = fit_type_bg

        self._num_pars = 0
        self._par_index_start_signal = 0
        self._par_index_start_bg = 0

        self._par_index_amplitude = 0
        self._par_index_mean = 0
        self._par_index_width = 0

        self._generate_fit_strings()
        self._generate_labels()

    def _generate_fit_strings(self):

        fit_par_index = 0

        self._par_index_start_signal = fit_par_index
        if self._fit_type_signal == FitTypeSignal.GAUS:
            self._fit_string_signal = f"[{fit_par_index}]*Numba::gaus(x, [{fit_par_index+1}], [{fit_par_index + 2}])"
            fit_par_index += 3

            self._par_index_amplitude = 0
            self._par_index_mean = 1
            self._par_index_width = 2

        elif self._fit_type_signal == FitTypeSignal.CRYSTAL:
            self._fit_string_signal = (
                f"[{fit_par_index}]*Numba::double_crystal(x, [{fit_par_index + 1}], [{fit_par_index + 2}], "
                f"[{fit_par_index + 3}], [{fit_par_index + 4}], [{fit_par_index + 5}], [{fit_par_index + 6}])"
            )

            # happens to be the same as the Gaussian case
            self._par_index_amplitude = 0
            self._par_index_mean = 1
            self._par_index_width = 2

            fit_par_index += 7

        # Background component
        self._par_index_start_bg = fit_par_index
        if self._fit_type_bg == FitTypeBg.POL1:
            self._fit_string_bg = f"pol1(0)"
            self._fit_string_bg_tmp = f"pol1({fit_par_index})"
            fit_par_index += 2
        elif self._fit_type_bg == FitTypeBg.POL2:
            self._fit_string_bg = f"pol2(0)"
            self._fit_string_bg_tmp = f"pol2({fit_par_index})"
            fit_par_index += 3

        # Total fit string
        self._fit_string_total = (
            f"{self._fit_string_signal} + {self._fit_string_bg_tmp}"
        )

        self._num_pars = fit_par_index

    def _generate_labels(self):
        self._fit_label_total = (
            f"{self._fit_type_signal.value} + {self._fit_type_bg.value}"
        )
        self._fit_label_signal = self._fit_type_signal.value
        self._fit_label_bg = self._fit_type_bg.value

    @property
    def fit_string_total(self):
        return self._fit_string_total

    @property
    def fit_string_signal(self):
        return self._fit_string_signal

    @property
    def fit_string_bg(self):
        return self._fit_string_bg

    @property
    def fit_label_total(self):
        return self._fit_label_total

    @property
    def fit_label_signal(self):
        return self._fit_label_signal

    @property
    def fit_label_bg(self):
        return self._fit_label_bg

    @property
    def num_pars(self) -> int:
        return self._num_pars

    @property
    def par_index_start_signal(self) -> int:
        return self._par_index_start_signal

    @property
    def par_index_start_bg(self) -> int:
        return self._par_index_start_bg

    @property
    def par_index_amplitude(self) -> int:
        return self._par_index_amplitude

    @property
    def par_index_mean(self) -> int:
        return self._par_index_mean

    @property
    def par_index_width(self) -> int:
        return self._par_index_width


class MassFit:
    def __init__(
        self,
        hist: rt.TH1D,
        fit_type_signal: FitTypeSignal,
        fit_type_bg: FitTypeBg,
        starting_pars: list[StartingPar],
        range_signal: RangeMassK0,
        range_sideband: RangeMassK0,
        fit_range: tuple[float, float] | None = None,
    ) -> None:

        self._id = uuid4()
        self._hist = hist.Clone(f"{self._id}_hist")
        self._fit_type_signal = fit_type_signal
        self._fit_type_bg = fit_type_bg
        self._starting_pars = starting_pars

        self._range_signal = range_signal
        self._range_sideband = range_sideband

        # Use histogram range if fit_range not provided
        if fit_range is None:
            self._fit_range = (
                self._hist.GetXaxis().GetXmin(),
                self._hist.GetXaxis().GetXmax(),
            )
        else:
            self._fit_range = fit_range

        self._fit_amplitude = 0
        self._fit_mean = 0
        self._fit_width = 0

        self._purity = 0
        self._missed_fraction = 0

        self._fit_descriptor = FitDescriptor(
            fit_type_signal=self._fit_type_signal, fit_type_bg=self._fit_type_bg
        )

        assert self._fit_descriptor.num_pars == len(self._starting_pars), (
            f"Number of starting parameters ({len(self._starting_pars)}) does not match "
            f"the number of parameters in the fit descriptor ({self._fit_descriptor.num_pars})."
        )

        self._fit_total = None
        self._fit_signal = None
        self._fit_bg = None

        self._configure_fits()
        self._fit_hist()
        self._compute_parameters()

    def _configure_fits(self):
        self._configure_fit_total()
        self._configure_fit_signal()
        self._configure_fit_bg()

    def _configure_fit_total(self):
        fit_total = rt.TF1(
            f"{self._id}_fit_total",
            self._fit_descriptor.fit_string_total,
            *self._fit_range,
        )

        for par_index, starting_par in enumerate(self._starting_pars):
            if starting_par.limits is not None:
                fit_total.SetParLimits(par_index, *starting_par.limits)
            if starting_par.fixed:
                fit_total.FixParameter(par_index, starting_par.value)
            else:
                fit_total.SetParameter(par_index, starting_par.value)

        self._fit_total = fit_total
        # for now we fix the BG parameters to straight line between fit range
        if self._fit_type_bg == FitTypeBg.POL1:
            self._overwrite_fit_parameters_bg()

        self._fit_total.SetNpx(1000)

    def _overwrite_fit_parameters_bg(self):
        bg_bin_low = self._hist.FindBin(self._fit_range[0])
        bg_bin_high = self._hist.FindBin(self._fit_range[1])

        bg_point_low = [self._fit_range[0], self._hist.GetBinContent(bg_bin_low)]
        bg_point_high = [self._fit_range[1], self._hist.GetBinContent(bg_bin_high)]

        bg_starting_params = get_straight_line(bg_point_low, bg_point_high)

        self._fit_total.FixParameter(
            self._fit_descriptor.par_index_start_bg,
            bg_starting_params[0],
        )
        self._fit_total.FixParameter(
            self._fit_descriptor.par_index_start_bg + 1,
            bg_starting_params[1],
        )

    def _configure_fit_signal(self):
        self._fit_signal = rt.TF1(
            f"{self._id}_fit_signal",
            self._fit_descriptor.fit_string_signal,
            *self._fit_range,
        )
        self._fit_signal.SetNpx(1000)

    def _configure_fit_bg(self):
        self._fit_bg = rt.TF1(
            f"{self._id}_fit_bg",
            self._fit_descriptor.fit_string_bg,
            *self._fit_range,
        )
        self._fit_bg.SetNpx(1000)

    def _fit_hist(self):
        self._fit_result = self._hist.Fit(
            self._fit_total,
            "RBS",
            "",
            *self._fit_range,
        )

        # Signal shape parameters
        for par_index in range(
            self._fit_descriptor.par_index_start_signal,
            self._fit_descriptor.par_index_start_bg,
        ):
            self._fit_signal.SetParameter(
                par_index - self._fit_descriptor.par_index_start_signal,
                self._fit_total.GetParameter(par_index),
            )

        # Update background component parameters
        for par_index in range(
            self._fit_descriptor.par_index_start_bg, self._fit_descriptor.num_pars
        ):
            self._fit_bg.SetParameter(
                par_index - self._fit_descriptor.par_index_start_bg,
                self._fit_total.GetParameter(par_index),
            )

    def _compute_parameters(self):
        lower_bound = (
            self.fit_mean.value + self.fit_width.value * self._range_signal.low
        )
        upper_bound = (
            self.fit_mean.value + self.fit_width.value * self._range_signal.high
        )
        self._purity = self._fit_signal.Integral(
            lower_bound, upper_bound
        ) / self._fit_total.Integral(lower_bound, upper_bound)

        self._missed_fraction = self._fit_signal.Integral(
            lower_bound, upper_bound
        ) / self._fit_signal.Integral(*self._fit_range)

    @property
    def fit_mean(self) -> Observable:
        """Extract the mean from the fit."""
        mean = self._fit_total.GetParameter(self._fit_descriptor.par_index_mean)
        mean_error = self._fit_total.GetParError(self._fit_descriptor.par_index_mean)
        return Observable(
            value=mean,
            stat_error=mean_error,
        )

    @property
    def fit_width(self) -> Observable:
        """Extract the width (sigma) from the fit."""

        width = self._fit_total.GetParameter(self._fit_descriptor.par_index_width)
        width_error = self._fit_total.GetParError(self._fit_descriptor.par_index_width)
        return Observable(
            value=width,
            stat_error=width_error,
        )

    @property
    def fit_amplitude(self) -> Observable:
        """Extract the amplitude from the fit."""

        amplitude = self._fit_total.GetParameter(
            self._fit_descriptor.par_index_amplitude
        )
        amplitude_error = self._fit_total.GetParError(
            self._fit_descriptor.par_index_amplitude
        )
        self._fit_amplitude = amplitude
        return Observable(
            value=amplitude,
            stat_error=amplitude_error,
        )

    @property
    def fit_total(self) -> rt.TF1:
        return self._fit_total

    @property
    def hist(self) -> rt.TH1D:
        return self._hist

    @property
    def fit_signal(self) -> rt.TF1:
        return self._fit_signal

    @property
    def fit_bg(self) -> rt.TF1:
        return self._fit_bg

    @property
    def fit_label_total(self) -> str:
        return self._fit_descriptor.fit_label_total

    @property
    def fit_label_signal(self) -> str:
        return self._fit_descriptor.fit_label_signal

    @property
    def fit_label_bg(self) -> str:
        return self._fit_descriptor.fit_label_bg

    @property
    def fit_result(self) -> rt.TFitResultPtr:
        return self._fit_result

    @property
    def chi2_ndf(self) -> float:
        """Return chi2/ndf of the fit."""
        return self._fit_total.GetChisquare() / self._fit_total.GetNDF()

    @property
    def range_signal(self) -> RangeMassK0:
        """Return the range of the signal region."""
        return self._range_signal

    @property
    def range_sideband(self) -> RangeMassK0:
        """Return the range of the sideband region."""
        return self._range_sideband

    @property
    def purity(self) -> float:
        """Return the purity of the signal region."""
        return self._purity

    @property
    def missed_fraction(self) -> float:
        return self._missed_fraction
