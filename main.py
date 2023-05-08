import matplotlib.pyplot as plt
from utils import *
# Analytical solution ##########################################################
_ksi = ksi(MASS, k, c)
_w_n = w_n(MASS, k)
_w_d = w_d(_ksi, _w_n)
consts = constants(X_0, X_DOT_0, _ksi, _w_n)


_x_h = x_h(t, _ksi, _w_n, _w_d, consts)
_x_p = x_p(t, F_0, MASS, k, c, w)

x = analytical_solution(_x_h, _x_p)
# Conv #########################################################################
# MDF ##########################################################################

# x_t_old = X_0

# x_t = X_DOT_0 * DELTA_T + x_t_old
# response = [x_t_old, x_t]
# for time in t[2:]:
#     A = k * x_t
#     B = c * ((x_t - x_t_old) / DELTA_T)

#     x_next = (
#         x_t / MASS * (2 * MASS - DELTA_T * c - DELTA_T**2 * k)
#         + x_t_old * (DELTA_T * c / MASS - 1)
#         - DELTA_T**2 / MASS * force(time)
#     )

#     response.append(x_next)

#     x_old = x_t
#     x_t = x_next

plt.plot(t, x.values())
plt.grid()

plt.show()
