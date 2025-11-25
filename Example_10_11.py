import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sv_from_coe import sv_from_coe
from lunar_position import lunar_position

# Global variable for Julian day
JD = None
mu = 398600
hours = 3600
days = 24 * hours
deg = np.pi / 180
mu3 = 4903
RE = 6378

def example_10_11():
    """
    This function solves Example 10.11 by using MATLAB’s ode45 to integrate
    Equations 10.84, the Gauss variational equations, for a lunar gravitational perturbation.
    """
    global JD

    # Conversion factors
    hours = 3600
    days = 24 * hours
    deg = np.pi / 180

    # Constants
    mu = 398600
    mu3 = 4903
    RE = 6378

    # Initial data for each of the three given orbits
    orbits = {
        'GEO': {'a0': 42164, 'e0': 0.0001, 'w0': 0, 'RA0': 0, 'i0': 1 * deg, 'TA0': 0, 'JD0': 2454283},
        'HEO': {'a0': 26553.4, 'e0': 0.741, 'w0': 270, 'RA0': 0, 'i0': 63.4 * deg, 'TA0': 0, 'JD0': 2454283},
        'LEO': {'a0': 6678.136, 'e0': 0.01, 'w0': 0, 'RA0': 0, 'i0': 28.5 * deg, 'TA0': 0, 'JD0': 2454283}
    }

    for n, data in orbits.items():
        global JD0
        JD0 = data['JD0']
        solveit(data, n, days, deg)

def solveit(data, n, days, deg):
    """
    Calculations and plots common to all of the orbits
    """
    # Initial orbital parameters (calculated from the given data)
    mu = 398600
    h0 = np.sqrt(mu * data['a0'] * (1 - data['e0'] ** 2))
    T0 = 2 * np.pi / np.sqrt(mu) * data['a0'] ** 1.5
    rp0 = h0 ** 2 / mu / (1 + data['e0'])
    ra0 = h0 ** 2 / mu / (1 - data['e0'])

    # Store initial orbital elements in the vector coe0
    coe0 = [h0, data['e0'], data['RA0'], data['i0'], data['w0'], data['TA0']]

    # Use solve_ivp to integrate the Gauss variational equations with lunar gravity as the perturbation
    t0 = 0
    tf = 60 * days
    y0 = coe0
    nout = 400
    tspan = np.linspace(t0, tf, nout)
    options = {'rtol': 1.e-8, 'atol': 1.e-8}

    sol = solve_ivp(rates, [t0, tf], y0, t_eval=tspan, method='RK45', **options)
    t = sol.t
    y = sol.y.T

    # Time histories of the right ascension, inclination, and argument of perigee
    RA = y[:, 2]
    i = y[:, 3]
    w = y[:, 4]

    # Smooth the data to eliminate short period variations
    RA = rsmooth(RA)
    i = rsmooth(i)
    w = rsmooth(w)

    plt.figure(n)
    plt.subplot(1, 3, 1)
    plt.plot(t / days, (RA - data['RA0']) / deg)
    plt.title('Right Ascension vs Time')
    plt.xlabel('t (days)')
    plt.ylabel('Ω (deg)')
    plt.axis('tight')

    plt.subplot(1, 3, 2)
    plt.plot(t / days, (i - data['i0']) / deg)
    plt.title('Inclination vs Time')
    plt.xlabel('t (days)')
    plt.ylabel('i (deg)')
    plt.axis('tight')

    plt.subplot(1, 3, 3)
    plt.plot(t / days, (w - data['w0']) / deg)
    plt.title('Argument of Perigee vs Time')
    plt.xlabel('t (days)')
    plt.ylabel('ω (deg)')
    plt.axis('tight')

    plt.show()

def rates(t, f):
    """
    Compute the rates for the ODE solver
    """
    global JD

    # Extract the orbital elements at time t
    h, e, RA, i, w, TA = f
    phi = w + TA  # argument of latitude

    # Obtain the state vector at time t from Algorithm 4.5
    coe = [h, e, RA, i, w, TA]
    R, V = sv_from_coe(coe, mu)

    # Obtain the unit vectors of the rsw system
    r = np.linalg.norm(R)
    ur = R / r
    H = np.cross(R, V)
    uh = H / np.linalg.norm(H)
    s = np.cross(uh, ur)
    us = s / np.linalg.norm(s)

    # Update the Julian Day
    JD = JD0 + t / days

    # Find and normalize the position vector of the moon
    R_m = lunar_position(JD)
    r_m = np.linalg.norm(R_m)
    R_rel = R_m - R
    r_rel = np.linalg.norm(R_rel)

    # See Appendix F
    q = np.dot(R, (2 * R_m - R)) / r_m ** 2
    F = (q ** 2 - 3 * q + 3) * q / (1 + (1 - q) ** 1.5)

    # Gravitational perturbation of the moon
    ap = mu3 / r_rel ** 3 * (F * R_m - R)

    # Perturbation components in the rsw system
    apr = np.dot(ap, ur)
    aps = np.dot(ap, us)
    aph = np.dot(ap, uh)

    # Gauss variational equations
    hdot = r * aps
    edot = h / mu * np.sin(TA) * apr + 1 / mu / h * ((h ** 2 + mu * r) * np.cos(TA) + mu * e * r) * aps
    RAdot = r / h / np.sin(i) * np.sin(phi) * aph
    idot = r / h * np.cos(phi) * aph
    wdot = -h * np.cos(TA) / mu / e * apr + (h ** 2 + mu * r) / mu / e / h * np.sin(TA) * aps - r * np.sin(phi) / h / np.tan(i) * aph
    TAdot = h / r ** 2 + 1 / e / h * (h ** 2 / mu * np.cos(TA) * apr - (r + h ** 2 / mu) * np.sin(TA) * aps)

    # Return rates to solve_ivp in the array dfdt
    dfdt = np.array([hdot, edot, RAdot, idot, wdot, TAdot])
    return dfdt

def rsmooth(y):
    """
    Function to smooth data to remove short period variations
    """
    # Implement a simple moving average smoothing as a placeholder
    window_size = 10
    y_smooth = np.convolve(y, np.ones(window_size) / window_size, mode='valid')
    return np.concatenate([y[:window_size-1], y_smooth])

example_10_11()
