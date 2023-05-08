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


def constants(x_0: float, x_dot_0: float, ksi: float, w_n: float) -> Tuple[float]:
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
    c_1 = x_0
    c_2 = (x_dot_0 + ksi * w_n * x_0) / (np.sqrt(1 - ksi ** 2) * w_n)
    return (c_1, c_2)


def x_h(t: List[float], ksi: float, w_n: float, w_d: float, constants: Tuple[int]) -> List[float]:
    c_1, c_2 = constants

    return [np.exp(-ksi * w_n * time) * (c_1 * np.cos(w_d * time) + c_2 * np.sin(w_d * time)) for time in t]


def x_p(t: List[float], F_0: float, mass: float, k: float, c: float, w: float) -> List[float]:

    C = F_0 / np.sqrt((k - mass * w ** 2) ** 2 + (c * w) ** 2)

    return C * np.cos(w * t - np.arctan(c * w / (k - mass * w ** 2)))


def analytical_solution(x_h: np.array, x_p: np.array) -> Dict[float, float]:
    resp = [x_h[i] + x_p[i] for i in range(len(x_h))]
    x = {t: resp[i] for i, t in enumerate(np.arange(0, t_f+DELTA_T, DELTA_T))}

    return x


def conv_solution(f: List[float], g: List[float]):

    return np.convolve(f, g)


def x_dot_t(x, t: List[float], delta_t: float) -> List[float]:

    return [(x(time) - x(time - delta_t)) / delta_t for time in t]


def x_ddot_t(x, t: List[float], delta_t: float) -> List[float]:

    return [(x(time + delta_t) - 2 * x(time) + x(time - delta_t)) / (delta_t ** 2) for time in t]
