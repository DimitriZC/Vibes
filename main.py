import matplotlib.pyplot as plt
from utils import *

# Analytical solution ##########################################################
_ksi = ksi(MASS, k, c)
_w_n = w_n(MASS, k)
_w_d = w_d(_ksi, _w_n)
_r = w/_w_n
consts = constants(X_0, X_DOT_0, _ksi, _w_n, MASS)

#solution using deltaT=0,02
_x_h = x_h(t, _ksi, _w_n, _w_d, consts)
Amp_xp, phi_xp, _x_p = x_p(t, F_0, MASS, k, c, w)  #Amp_xp refers to amplitude of the mass caused by the force
analyt_sol, dict_sol = analytical_solution(_x_h, _x_p, deltaT=DELTA_T)

#solution using deltaT=0,002
_x_h2 = x_h(t2, _ksi, _w_n, _w_d, consts)
Amp_xp2, phi_xp2, _x_p2 = x_p(t2, F_0, MASS, k, c, w)  #Amp_xp refers to amplitude of the mass caused by the force
analyt_sol2, dict_sol2 = analytical_solution(_x_h2, _x_p2, deltaT=DELTA_T2)


# Solution by convolution integral #############################################
conv_t = np.array(t)
conv_t2 = np.array(t2)
conv_solut = conv_solution(conv_t, _ksi, _w_n, _w_d, consts, deltaT=DELTA_T)  #for deltaT=0,02
conv_solut2 = conv_solution(conv_t2, _ksi, _w_n, _w_d, consts, deltaT=DELTA_T2)  #for deltaT=0,002

# MDF ##########################################################################

FDM_solution, time_array1 = FDM_solver(deltat = DELTA_T, t_final=t_f)
FDM_solution2, time_array2 = FDM_solver(deltat = DELTA_T2, t_final=t_f)

# print([(round(conv_solut[i], 5), round(FDM_solution[i], 5), round(analyt_sol[i], 5)) for i in range(10)]) #comparing first 10 points from convolution, FDM and analytical methods


# Plots ########################################################################
def plot_3_methods_comparison(time, anal_solution, convol_solution, fdifferences_solution):
    plt.figure(figsize=(10,4))
    plt.ylabel('Posição [m]')
    plt.xlabel('Tempo [s]')
    plt.plot(time, anal_solution, label = 'Solução analítica', color = 'orange') #plots FDM solution
    plt.plot(time, convol_solution, label = 'Integral de convolução', color = 'g')  #plots convolution solution
    plt.plot(time, fdifferences_solution, label = 'Diferenças finitas', color = 'b') #plots FDM solution
    plt.legend()
    plt.grid()
    plt.show()
    return 

def plot_deltaT_diff(solu_deltat1, time_array1, solu_deltat2, time_array2): 
    plt.figure(figsize=(10,4))
    plt.ylabel('Posição [m]')
    plt.xlabel('Tempo [s]')
    plt.plot(time_array1, solu_deltat1, label = r'$\Delta t = 0,02s$', color = 'orange') #plots FDM solution
    plt.plot(time_array2, solu_deltat2, label = r'$\Delta t = 0,002s$', color = 'b')  #plots convolution solution
    plt.legend()
    plt.grid()
    plt.show()
    return 

def plot_FDM_error(FDM_solution, time_array1, exact_solution, time_array2):
  
    FDM_error = np.array(FDM_solution) - np.array(exact_solution)
    FDM_approx_deriv = np.diff(FDM_solution)/(time_array1[1]-time_array1[0])
    fig, ax1 = plt.subplots(figsize=(10,4))

    ax1.set_ylabel('Posição [m]')
    ax1.set_xlabel('Tempo [s]')
    ax1.plot(time_array1, FDM_error, label = 'Erro de diferenças finitas', color = 'orange') #plots FDM solution
    ax1.legend(loc='upper left')
    ax1.tick_params('y', colors='orange')

    ax2 = ax1.twinx()
    ax2.plot(time_array2[:-1], FDM_approx_deriv, label = 'Derivada de diferenças finitas', color = 'b')  #plots convolution solution
    ax2.legend(loc='upper right')
    ax2.tick_params('y', colors='b')

    plt.grid()
    plt.show()
    return 

#### uncomment to plot three methods for deltaT=0,02
plot_3_methods_comparison(t[:-1], analyt_sol[:-1], conv_solut[:len(t)-1], FDM_solution[:len(t)-1])

#### uncomment to plot three methods for deltaT=0,002
plot_3_methods_comparison(t2[:-1], analyt_sol2[:-1], conv_solut2[:len(t2)-1], FDM_solution2[:len(t2)-1])

#### uncomment to plot analytical method for different deltaT's
plot_deltaT_diff(analyt_sol[:len(t)-1], t[:-1], analyt_sol2[:len(t2)-1], t2[:-1])

#### uncomment to plot FDM method for different deltaT's
plot_deltaT_diff(FDM_solution[:len(time_array1)-1], time_array1[:-1], FDM_solution2[:len(time_array2)-1], time_array2[:-1])

#### uncomment to plot convolution method for different deltaT's
plot_deltaT_diff(conv_solut[:len(conv_t)-1], conv_t[:-1], conv_solut2[:len(conv_t2)-1], conv_t2[:-1])

#### uncomment to plot the error from the FDM method and its derivatives
# plot_FDM_error(FDM_solution, t, analyt_sol, t)
# print(maximum_diff(t, analyt_sol, FDM_solution))   #prints the maximum difference found in two different methods