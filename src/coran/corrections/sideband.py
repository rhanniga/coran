from uuid import uuid4, UUID

import ROOT as rt


class SidebandCorrector:
    def __init__(
        self,
        signal_2d: rt.TH2D,
        sideband_2d: rt.TH2D,
        purity: float,
        missed_fraction: float,
        branching_ratio: float = 0.692,
    ):
        """This class corrects for the combinatorial background using the sideband technique,
        as well as a few smaller corrections that are too trivial to isolate (missed fraction
        from signal region and the branching ratio of the V0)"""

        self._id = uuid4()
        self._signal_2d = signal_2d.Clone(f"{self._id}_signal_2d")
        self._sideband_2d = sideband_2d.Clone(f"{self._id}_sideband_2d")

        self._corrected_2d = signal_2d.Clone(f"{self._id}_corrected_2d")

        self._purity = purity
        self._missed_fraction = missed_fraction
        self._branching_ratio = branching_ratio

        self._process_corrections()

    def _process_corrections(self):
        self._sideband_2d.Scale(1 / self._sideband_2d.Integral())

        bg_integral = (1 - self._purity) * self._corrected_2d.Integral()

        self._corrected_2d.Add(self._sideband_2d, -bg_integral)
        self._corrected_2d.Scale(1 / self._missed_fraction)
        self._corrected_2d.Scale(1 / self._branching_ratio)

    @property
    def sideband_2d(self) -> rt.TH2D:
        return self._sideband_2d

    @property
    def signal_2d(self) -> rt.TH2D:
        return self._signal_2d

    @property
    def corrected_2d(self) -> rt.TH2D:
        return self._corrected_2d
