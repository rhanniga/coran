EPSILON = 1e-6

INPUT_FILE = "input/lhc16q_final.root"
INPUT_LIST = "h-k0"

BRANCHING_RATIO_K0 = 0.6920

NUM_BINS_Z_VERTEX = 10

RANGE_PT_TRIG = (4.0, 8.0 - EPSILON)
RANGE_DELTA_ETA = (-1.2, 1.2 - EPSILON)
RANGE_FIT_MASS = (0.44, 0.56 - EPSILON)

RANGES_MULTIPLICITY = {
    "mult_0_20" :  (0.0, 20.0 - EPSILON),
    "mult_20_50" :(20.0, 50.0 - EPSILON),
    "mult_50_80" : (50.0, 80.0 - EPSILON),
}

RANGES_PT_ASSOCIATED = {
    "pt_10_15": (1.0, 1.5 - EPSILON),
    "pt_15_20": (1.5, 2.0 - EPSILON),
    "pt_20_25": (2.0, 2.5 - EPSILON),
    "pt_25_30": (2.5, 3.0 - EPSILON),
    "pt_30_40": (3.0, 4.0 - EPSILON),
    "pt_15_25": (1.5, 2.5 - EPSILON),
    "pt_25_40": (2.5, 4.0 - EPSILON),
}

RANGE_SIGNAL = (-5.0, 5.0 - EPSILON)

VARIATIONS_SIGNAL = {
    "wide": (-8.0, 8.0 - EPSILON),
    "wider": (-10.0, 10.0 - EPSILON),
    "narrow":(-3.0, 3.0 - EPSILON),
}

RANGE_SIDEBAND = (-16.5, -6.5 - EPSILON)

VARIATIONS_SIDEBAND = {
    "shifted_left": (-19.5, -9.5 - EPSILON),
    "shifted_right": (6.5, 16.5 - EPSILON),
    "shifted_more_right": (9.5, 19.5 - EPSILON),
}



AXIS_HADRON_PT = 0
AXIS_HADRON_PHI = 1
AXIS_HADRON_ETA = 2
AXIS_HADRON_ZVTX = 3
AXIS_HADRON_MULT = 4

# For single particle v0 distributions
AXIS_V0_PT = 0
AXIS_V0_PHI = 1
AXIS_V0_ETA = 2
AXIS_V0_MASS = 3
AXIS_V0_MULT = 4

# For h-h correlation distributions
AXIS_H_H_PT_TRIGGER = 0
AXIS_H_H_PT_ASSOCIATED = 1
AXIS_H_H_DELTA_PHI = 2
AXIS_H_H_DELTA_ETA = 3
AXIS_H_H_ZVTX = 4
AXIS_H_H_MULT = 5

# For h-k correlation distributions
AXIS_H_K_PT_TRIGGER = 0
AXIS_H_K_PT_ASSOCIATED = 1
AXIS_H_K_DELTA_PHI = 2
AXIS_H_K_DELTA_ETA = 3
AXIS_H_K_MASS = 4
AXIS_H_K_ZVTX = 5
AXIS_H_K_MULT = 6
