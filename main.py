import matplotlib.pyplot as plt
from utils import *

# Analytical solution ##########################################################
_ksi = ksi(MASS, k, c)
_w_n = w_n(MASS, k)
_w_d = w_d(_ksi, _w_n)
_r = w/_w_n
consts = constants(X_0, X_DOT_0, _ksi, _w_n, MASS)

_x_h = x_h(t, _ksi, _w_n, _w_d, consts)
Amp_xp, phi_xp, _x_p = x_p(t, F_0, MASS, k, c, w)  #Amp_xp refers to amplitude of the mass caused by the force
plot_x, dict_x = analytical_solution(_x_h, _x_p)
# print('resposta em x({0}): {1}'.format(0, dict_x[0]))


# plt.plot(t[:-1], _x_p[:-1])
# plt.grid()
# plt.show()

# Conv #########################################################################

conv_t = np.array(t)
conv_xp = conv_solution(conv_t, _ksi, _w_n, _w_d, consts)  

print([(conv_xp[i].round(5), plot_x[i].round(5)) for i in range(10)]) #comparing first 10 points from convolution method and analytical
plt.plot(t[:-1], conv_xp[:len(t)-1])
plt.grid()

plt.show()

# MDF ##########################################################################

# x_t_old = X_0

# x_t = X_DOT_0 * DELTA_T + x_t_old
# response = [x_t_old, x_t]
# for time in t[2:]:
#     A = k * x_t
#     B = c * ((x_t - x_t_old) / DELTA_T)

#     x_next = (
#         (x_t / MASS) * (2 * MASS - DELTA_T * c - DELTA_T**2 * k)
#         + x_t_old * (DELTA_T * c / MASS - 1)
#         - DELTA_T**2 / MASS * force(time)
#     )

#     response.append(x_next)

#     x_t_old = x_t
#     x_t = x_next

# plt.plot(t[:-1], response[:-1])
# plt.grid()

# plt.show()

"""
- [x] Plotar solução Homogenea;
- [x] Plotar Solução Particular;
- [] Revisar Solução analítica(constantes);
- [] Colocar DeltaT diferente para solução analítica
- [] Melhorar os plots;
- [] Entender a convolução;
=>
- [] Escrever o relatório;
"""
