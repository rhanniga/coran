EPSILON = 1e-6

INPUT_FILENAME = "lhc16q_final.root"
INPUT_LIST = "h-k0"

RANGES_MULTIPLICITY = [
    (0.0, 20.0 - EPSILON),
    (20.0, 50.0 - EPSILON),
    (50.0, 80.0 - EPSILON),
]

RANGES_PT_ASSOCIATED = [
    (1.0, 1.5 - EPSILON),
    (1.5, 2.0 - EPSILON),
    (2.0, 2.5 - EPSILON),
    (2.5, 3.0 - EPSILON),
    (3.0, 4.0 - EPSILON),
]

# For single particle hadron distributions
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
