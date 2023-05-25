import numpy as np

k = 1_000  # [N / m]
c = 5  # [N * s / m]
MASS = 1  # [kg]

X_0 = 0.7  # [m]
X_DOT_0 = 30  # [m / s]

DELTA_T = 0.002 # [s]
t_f = 5 # [s]
F_0 = -100
w = 50

t = [i for i in np.arange(0, 5+DELTA_T, DELTA_T)]
