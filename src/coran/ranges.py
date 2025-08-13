from abc import ABC, abstractmethod


class RangeBase(ABC):
    def __init__(self, low: float, high: float):
        self._low = low
        self._high = high

    @abstractmethod
    def __repr__(self) -> str:
        pass

    @property
    def low(self) -> float:
        return self._low

    @property
    def high(self) -> float:
        return self._high


class RangeMult(RangeBase):
    def __repr__(self) -> str:
        return f"{self._low:.1f}#minus{self._high:.1f}% V0A"


class RangePtAssoc(RangeBase):
    def __repr__(self) -> str:
        return f"{self._low:.1f} < #it{{p}}_{{T, assoc}} {self._high:.1f} GeV/#it{{c}}"


class RangePtTrig(RangeBase):
    def __repr__(self) -> str:
        return f"{self._low:.1f} < #it{{p}}_{{T, trig}} {self._high:.1f} GeV/#it{{c}}"


class RangeDeltaEta(RangeBase):
    def __repr__(self) -> str:
        if abs(self._low) == abs(self._high):
            return f"|#Delta#it{{#eta}}| < {abs(self._low):.1f}"
        return f"{self._low:.1f} < #Delta#it{{#eta}} < {self._high:.1f}"


class RangeMassK0(RangeBase):
    def __repr__(self) -> str:
        return f"{self._low:.3f} < #it{{M}}_{{#pi^{{+}}#pi^{{-}}}} < {self._high:.3f} GeV/#it{{c}}^{2}"
