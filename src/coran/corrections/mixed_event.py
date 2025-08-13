from enum import Enum
from uuid import uuid4, UUID

import ROOT as rt


class ScaleType(Enum):
    ZERO = "#Delta#varphi, #Delta#eta = 0, 0"
    PI = "#Delta#varphi, #Delta#eta = #pi, 0"
    AVG = "Avg. #Delta#varphi# along #Delta#eta = 0"


class AcceptanceCorrector:
    def __init__(
        self,
        same: rt.TH3D,
        mixed: rt.TH3D,
        scale_type: ScaleType = ScaleType.ZERO,
        axis_delta_eta: int = 1,
        axis_delta_phi: int = 2,
        axis_z_vertex: int = 3,
        num_bins_z_vertex: int = 10,
    ):

        self._id = uuid4()
        self._same = same.Clone(f"{self._id}_same")
        self._mixed = mixed.Clone(f"{self._id}_mixed")

        self._same_2d = same.Project3D("xye")
        self._mixed_2d = mixed.Project3D("xye")
        self._corrected_2d = None

        self._scale_type = scale_type
        self._axis_delta_eta = axis_delta_eta
        self._axis_delta_phi = axis_delta_phi
        self._axis_z_vertex = axis_z_vertex
        self._num_bins_z_vertex = num_bins_z_vertex

        self._do_correction()

    def _do_correction(self):
        for z_vertex_bin in range(1, self._num_bins_z_vertex + 1):
            self._same.GetZaxis().SetRange(z_vertex_bin, z_vertex_bin)
            same_2d_partial = self._same.Project3D("xye")

            self._mixed.GetZaxis().SetRange(z_vertex_bin, z_vertex_bin)
            mixed_2d_partial = self._mixed.Project3D("xye")

            if self._scale_type == ScaleType.ZERO:
                scale = 0.5 * (
                    mixed_2d_partial.Integral(
                        mixed_2d_partial.GetXaxis().FindBin(-0.01),  # xmin
                        mixed_2d_partial.GetXaxis().FindBin(0.01),  # xmax
                        mixed_2d_partial.GetYaxis().FindBin(0.0),  # ymin
                        mixed_2d_partial.GetYaxis().FindBin(0.0),
                    )
                )
            else:
                raise NotImplementedError(
                    f"Scale type {self._scale_type} is not implemented."
                )

            same_2d_partial.Divide(mixed_2d_partial)
            same_2d_partial.Scale(scale)

            if z_vertex_bin == 1:
                self._corrected_2d = same_2d_partial.Clone(f"{self._id}_same_2d")
            else:
                self._corrected_2d.Add(same_2d_partial)

    @property
    def mixed_2d(self) -> rt.TH2D:
        return self._mixed_2d

    @property
    def same_2d(self) -> rt.TH2D:
        return self._same_2d

    @property
    def corrected_2d(self) -> rt.TH2D:
        return self._corrected_2d
