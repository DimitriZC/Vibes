import numpy as np
from typing import Tuple, List, Dict
from .settings import *


def force(t: float) -> float:
    """This functions calculates the force in the time t

    Args:
        t (float): time

    Returns:
        float: Forece [N]
    """
    return F_0 * np.cos(w * t)


def ksi(mass: float, k: float, c: float) -> float:
    """This function calculates the damping ratio of the system

    Args:
        mass (float): Mass of the system [kg]
        k (float): spring constant system [N / m]
        c (float): damping coefficient [N * s / m]

    Returns:
        float: damping ratio [-]
    """
    return c / (2 * np.sqrt(k * mass))




def w_n(mass: float, k: float) -> float:
    """This functions calculates the natural frequency of the system

    Args:
        mass (float): Mass of the system [kg]
        k (float): spring constant of the system [N / m]

    Returns:
        float: natural frequency [rad / s]
    """
    return np.sqrt(k / mass)


def w_d(ksi: float, w_n):
    """This function calculates the damped natural frequency

    Args:
        ksi (float): damping ratio [-]
        w_n (_type_): natural frequency of the system [rad / s]

    Returns:
        float: damped natural frequency [rad / s]
    """
    return np.sqrt(1 - ksi**2) * w_n


def constants(x_0: float, x_dot_0: float, ksi: float, w_n: float, mass: float) -> Tuple[float]:
    """This function determinates the coefficients of the homogeneous solution
    of the motion equation

    Args:
        x_0 (float): initial position [m]
        x_dot_0 (float): initial velocity [m / s]
        ksi (float): damping ratio [-]
        w_n (float): natural frequency [rad / s]

    Returns:
        Tuple[float]: coefficients of the homogeneous solution
    """
    Amp = F_0 / (np.sqrt((k - mass * w ** 2) ** 2 + (c * w) ** 2))
    phi = np.arctan(c * w / (k - mass * w ** 2)) + np.pi
    c_1 = x_0 - Amp*np.cos(-phi)
    c_2 = (x_dot_0  + Amp* w * np.sin(-phi) + ksi * w_n * c_1) / (np.sqrt(1 - ksi ** 2) * w_n)
    return (c_1, c_2)


def x_h(t: List[float], ksi: float, w_n: float, w_d: float, constants: Tuple[int]) -> List[float]:
    c_1, c_2 = constants
    Amp_xh = np.sqrt(X_0 ** 2 + (ksi * w_n * X_0 / w_d)**2)
    phi_xh = np.arctan((X_DOT_0 + ksi * w_n * X_0) / (w_d * X_0))
    # return [Amp_xh * np.exp(-ksi * w_n * time) * np.cos(w_d*time - phi_xh) for time in t]   #desse jeito está dando errado
    return [np.exp(-ksi * w_n * time) * (c_1 * np.cos(w_d * time) + c_2 * np.sin(w_d * time)) for time in t]


def x_p(t: List[float], F_0: float, mass: float, k: float, c: float, w: float) -> List[float]:

    Amp = F_0 / (np.sqrt((k - mass * w ** 2) ** 2 + (c * w) ** 2))
    phi = np.arctan(c * w / (k - mass * w ** 2)) + np.pi

    return Amp, phi, [Amp * np.cos(w * time - phi) for time in t]



def analytical_solution(x_h: np.array, x_p: np.array, deltaT) -> Dict[float, float]:
    resp = [x_h[i] + x_p[i] for i in range(len(x_h))]
    x = {t: resp[i] for i, t in enumerate(np.arange(0, t_f+deltaT, deltaT))}

    return resp, x 



def conv_solution(tt: List[float], ksi: float, w_n: float, w_d: float, constants: Tuple[int], deltaT) -> List[float]:

    f_t = [F_0*np.cos(w * time) for time in tt] 

    g_t = [(np.exp(-ksi*w_n*time) * np.sin(w_d*time))/(MASS*w_d) for time in tt]
    xp_conv = np.convolve(f_t, g_t) * deltaT

    # to find C_1 and C_2 we do x_tot(0) = 0.7 and x_dot_tot(0) = 30, in equations:
    # x(0) = x_hom(0) + xp_conv(0) = 0.7
    # x_dot(0) = x_dot_hom(0) + xp_dot_conv(0) = 30, and we can approximate xp_dot_conv by: xp_dot_conv(x) = (xp_conv(t+dt) - xp_conv(t))/dt

    c_1_conv = X_0 - xp_conv[0]
    x_dot_0_xp_conv = (xp_conv[1] - xp_conv[0])/deltaT
    c_2_conv = (X_DOT_0 + c_1_conv * ksi * w_n + x_dot_0_xp_conv) / w_d   ####CHECK THIS EQUATION
    x_h_conv = [np.exp(-ksi * w_n * time) * (c_1_conv * np.cos(w_d * time) + c_2_conv * np.sin(w_d * time)) for time in tt]
    x_tot = np.add(xp_conv[:len(x_h_conv)], x_h_conv)
    # print('constants for convolution method:({0}, {1})'.format(c_1_conv, c_2_conv))
   


    return x_tot

def FDM_solver(deltat, t_final):
    tt = [i for i in np.arange(0, t_final+deltat, deltat)]
    x_t_old = X_0

    x_t = X_DOT_0 * deltat + x_t_old
    FDM_solution = [x_t_old, x_t]

    for time in tt[2:]:  
        kt = k * x_t
        ct = c * ((x_t - x_t_old) / deltat)

        x_next = (
    (force(time) - kt - ct) * (deltat**2 / MASS) 
        + 2 * x_t - x_t_old
        )

        FDM_solution.append(x_next)

        x_t_old = x_t
        x_t = x_next
    return FDM_solution, tt

def x_dot_t(x, t: List[float], delta_t: float) -> List[float]:

    return [(x(time) - x(time - delta_t)) / delta_t for time in t]


def x_ddot_t(x, t: List[float], delta_t: float) -> List[float]:

    return [(x(time + delta_t) - 2 * x(time) + x(time - delta_t)) / (delta_t ** 2) for time in t]


def maximum_diff(time, method1, method2):
    max_diff = 0
    temp_time = 0
    for i in range(len(time)):
        temp_diff = abs(method1[i] - method2[i])
        if temp_diff >= max_diff:
            max_diff = temp_diff
            temp_time = time[i]
    return temp_time, max_diff
