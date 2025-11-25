import numpy as np
from scipy.integrate import solve_ivp
from math import pi, sqrt
from coe_from_sv import coe_from_sv
from rv_from_r0v0_ta import rv_from_r0v0_ta
from rkf45 import rkf45

# Constants
mu = 398600  # Gravitational parameter (km^3/s^2)
RE = 6378  # Earth's radius (km)
g0 = 9.807  # Sea level acceleration of gravity (m/s^2)
deg = pi / 180  # Degree to radian conversion

# Input data
r0 = np.array([RE + 480, 0, 0])  # Initial position vector (km)
v0 = np.array([0, 7.7102, 0])  # Initial velocity vector (km/s)
t0 = 0  # Initial time (s)
t_burn = 261.1127  # Burn time (s)

m0 = 2000  # Initial spacecraft mass (kg)
T = 10  # Rated thrust of rocket engine (kN)
Isp = 300  # Specific impulse of rocket engine (s)

# Initial state vector
y0 = np.concatenate((r0, v0, [m0]))

# Function to compute acceleration vector
def rates(t, f):
    x = f[0]
    y = f[1]
    z = f[2]
    vx = f[3]
    vy = f[4]
    vz = f[5]
    m = f[6]

    r = np.linalg.norm([x, y, z])
    v = np.linalg.norm([vx, vy, vz])
    
    ax = -mu * x / r**3 + T / m * vx / v
    ay = -mu * y / r**3 + T / m * vy / v
    az = -mu * z / r**3 + T / m * vz / v  # Include gravity acceleration here
    mdot = -T * 1000 / g0 / Isp  # Rate of change of mass

    dfdt = np.array([vx, vy, vz, ax, ay, az, mdot])
    # transform dfdt
    return dfdt

# Numerical integration using rkf45
t_span = [t0, t_burn]
t, y = rkf45(rates, t_span, y0, 1e-16)

# Extracting results after burn
r1 = y[-1, :3]
v1 = y[-1, 3:6]
m1 = y[-1, 6]

# Compute orbital elements of the post-burn trajectory
coe = coe_from_sv(r1, v1, mu)

# Find state vector at apogee
e = coe[1]  # eccentricity
TA = coe[5]  # true anomaly (radians)
a = coe[6]  # semimajor axis (km)

if TA <= pi:
    dtheta = pi - TA
else:
    dtheta = 3 * pi - TA

ra, va = rv_from_r0v0_ta(r1, v1, dtheta/deg, mu)
rmax = np.linalg.norm(ra)

# Function to output results
def output(r0, v0, m0, r1, v1, m1, coe, ra, va, rmax):
    print('\n\n–––––––––––––––––––––––––––––––––––––––––––––––––––––')
    print('\nBefore ignition:')
    print(f' Mass = {m0} kg')
    print(' State vector:')
    print(f' r = [{r0[0]:.10g}, {r0[1]:.10g}, {r0[2]:.10g}] (km)')
    print(f' Radius = {np.linalg.norm(r0)} km')
    print(f' v = [{v0[0]:.10g}, {v0[1]:.10g}, {v0[2]:.10g}] (km/s)')
    print(f' Speed = {np.linalg.norm(v0)} km/s')
    print(f'\nThrust: {T} kN')
    print(f'Burn time: {t_burn:.6f} s')
    print(f'Mass after burn = {m1:.6E} kg\n')
    print('End-of-burn-state vector:')
    print(f' r = [{r1[0]:.10g}, {r1[1]:.10g}, {r1[2]:.10g}] (km)')
    print(f' Radius = {np.linalg.norm(r1)} km')
    print(f' v = [{v1[0]:.10g}, {v1[1]:.10g}, {v1[2]:.10g}] (km/s)')
    print(f' Speed = {np.linalg.norm(v1)} km/s\n')
    print('Post-burn trajectory:')
    print(f' Eccentricity = {coe[1]}')
    print(f' Semimajor axis = {coe[6]} km')
    print(' Apogee state vector:')
    print(f' r = [{ra[0]:.10E}, {ra[1]:.10E}, {ra[2]:.10E}] (km)')
    print(f' Radius = {np.linalg.norm(ra)} km')
    print(f' v = [{va[0]:.10E}, {va[1]:.10E}, {va[2]:.10E}] (km/s)')
    print(f' Speed = {np.linalg.norm(va)}')
    print('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')

# Output results
output(r0, v0, m0, r1, v1, m1, coe, ra, va, rmax)
