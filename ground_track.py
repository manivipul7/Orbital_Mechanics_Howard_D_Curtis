import numpy as np
import matplotlib.pyplot as plt
from kepler_E import kepler_E
from ra_and_dec_from_r import ra_and_dec_from_r
from sv_from_coe import sv_from_coe

# Global variables for storing right ascensions and declinations
ra = []
dec = []
n_curves = 0
RA = []
Dec = []

# Constants
deg = np.pi / 180
mu = 398600
J2 = 0.00108263
Re = 6378
we = (2 * np.pi + 2 * np.pi / 365.26) / (24 * 3600)

def ground_track():
    """
    This function computes the ground track of an Earth satellite based on specified orbital elements.
    It then plots the ground track and prints orbital data.
    """
    global ra, dec, n_curves, RA, Dec
    
    # Clear previous data
    ra = []
    dec = []
    RA = []
    Dec = []
    
    # Constants and data declaration for Example 4.12
    rP = 6700
    rA = 10000
    TAo = 230 * deg
    Wo = 270 * deg
    incl = 60 * deg
    wpo = 45 * deg
    n_periods = 3.25
    
    # Compute derived quantities
    a = (rA + rP) / 2
    T = 2 * np.pi / np.sqrt(mu) * a ** (3 / 2)
    e = (rA - rP) / (rA + rP)
    h = np.sqrt(mu * a * (1 - e ** 2))
    Eo = 2 * np.arctan(np.tan(TAo / 2) * np.sqrt((1 - e) / (1 + e)))
    Mo = Eo - e * np.sin(Eo)
    to = Mo * (T / (2 * np.pi))
    tf = to + n_periods * T
    fac = -3 / 2 * np.sqrt(mu) * J2 * Re ** 2 / ((1 - e ** 2) ** 2 * a ** (7 / 2))
    Wdot = fac * np.cos(incl)
    wpdot = fac * (5 / 2 * np.sin(incl) ** 2 - 2)
    
    # Find RA and Dec over the specified time interval
    find_ra_and_dec(to, tf, T, e, h, Mo, Wo, wpo, Wdot, wpdot, incl)
    
    # Form separate curves for plotting
    form_separate_curves()
    
    # Plot the ground track
    plot_ground_track()
    
    # Print orbital data
    print_orbital_data(h, e, a, rP, rA, T, incl, TAo, to, Wdot, Wo, wpdot, wpo)

def find_ra_and_dec(to, tf, T, e, h, Mo, Wo, wpo, Wdot, wpdot, incl):
    """
    Propagates the orbit over the specified time interval, transforming
    the position vector into the earth-fixed frame and computes the
    right ascension and declination histories.
    """
    global ra, dec
    
    times = np.linspace(to, tf, 1000)
    ra = []
    dec = []
    theta = 0
    
    for t in times:
        M = 2 * np.pi / T * t
        E = kepler_E(e, M)
        TA = 2 * np.arctan(np.tan(E / 2) * np.sqrt((1 + e) / (1 - e)))
        
        r = h ** 2 / mu / (1 + e * np.cos(TA)) * np.array([np.cos(TA), np.sin(TA), 0])
        
        W = Wo + Wdot * t
        wp = wpo + wpdot * t
        
        R1 = np.array([[np.cos(W), np.sin(W), 0],
                       [-np.sin(W), np.cos(W), 0],
                       [0, 0, 1]])
        
        R2 = np.array([[1, 0, 0],
                       [0, np.cos(incl), np.sin(incl)],
                       [0, -np.sin(incl), np.cos(incl)]])
        
        R3 = np.array([[np.cos(wp), np.sin(wp), 0],
                       [-np.sin(wp), np.cos(wp), 0],
                       [0, 0, 1]])
        
        QxX = np.dot(R3, np.dot(R2, R1))
        R = np.dot(QxX, r)
        
        theta = we * (t - to)
        Q = np.array([[np.cos(theta), np.sin(theta), 0],
                      [-np.sin(theta), np.cos(theta), 0],
                      [0, 0, 1]])
        
        r_rel = np.dot(Q, R)
        alpha, delta = ra_and_dec_from_r(r_rel)
        
        ra.append(alpha)
        dec.append(delta)

def form_separate_curves():
    """
    Breaks the ground track up into separate curves which start
    and terminate at right ascensions in the range [0, 360 deg].
    """
    global ra, dec, n_curves, RA, Dec
    
    tol = 100
    curve_no = 1
    n_curves = 1
    k = 0
    ra_prev = ra[0]
    
    RA.append([])
    Dec.append([])
    
    for i in range(len(ra)):
        if np.abs(ra[i] - ra_prev) > tol:
            curve_no += 1
            n_curves += 1
            k = 0
        
        k += 1
        RA.append([])
        Dec.append([])
        RA[curve_no].append(ra[i])
        Dec[curve_no].append(dec[i])
        ra_prev = ra[i]

def plot_ground_track():
    """
    Plots the ground track using Matplotlib.
    """
    global RA, Dec, n_curves
    
    plt.figure()
    plt.xlabel('East longitude (degrees)')
    plt.ylabel('Latitude (degrees)')
    plt.axis('equal')
    plt.grid(True)
    
    for i in range(1, n_curves + 1):
        plt.plot(RA[i], Dec[i])
    
    plt.axis([0, 360, -90, 90])
    plt.text(ra[0], dec[0], 'Start')
    plt.text(ra[-1], dec[-1], 'Finish')
    plt.axhline(0, color='black')
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    plt.show()

def print_orbital_data(h, e, a, rP, rA, T, incl, TAo, to, Wdot, Wo, wpdot, wpo):
    """
    Prints the orbital data.
    """
    coe = [h, e, Wo, incl, wpo, TAo]
    ro, vo = sv_from_coe(coe, mu)
    
    print('\n---------------------------------------------')
    print('\n Orbital Parameters:')
    print(f'\n Angular momentum (h): {h} km^2/s')
    print(f'\n Eccentricity (e): {e}')
    print(f'\n Semimajor axis (a): {a} km')
    print(f'\n Perigee radius (rP): {rP} km')
    print(f'\n Apogee radius (rA): {rA} km')
    print(f'\n Period (T): {T / 3600} hours')
    print(f'\n Inclination (incl): {np.rad2deg(incl)} deg')
    print(f'\n Initial true anomaly (TAo): {np.rad2deg(TAo)} deg')
    print(f'\n Time since perigee (to): {to / 3600} hours')
    print(f'\n Initial RA (Wo): {np.rad2deg(Wo)} deg')
    print(f'\n RA dot (Wdot): {np.rad2deg(Wdot) * T} deg/period')
    print(f'\n Initial argument of perigee (wpo): {np.rad2deg(wpo)} deg')
    print(f'\n Argument of perigee dot (wpdot): {np.rad2deg(wpdot) * T} deg/period')
    print(f'\n')
    print(f'\n Initial position vector (r0): {ro} km')
    print(f'\n Magnitude: {np.linalg.norm(ro)} km')
    print(f'\n Initial velocity vector (v0): {vo} km/s')
    print(f'\n Magnitude: {np.linalg.norm(vo)} km/s')
    print('\n---------------------------------------------')

ground_track()