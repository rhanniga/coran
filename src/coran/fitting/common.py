from dataclasses import dataclass


@dataclass
class StartingPar:
    value: float
    limits: tuple[float, float] | None = None
    fixed: bool = False


@dataclass
class Observable:
    value: float
    stat_error: float
    sys_error: float | None = None
