import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sv_from_coe import sv_from_coe
from atmosphere import atmosphere

# Conversion factors and constants
hours = 3600
days = 24 * hours
deg = np.pi / 180

mu = 398600  # Gravitational parameter (km^3/s^2)
RE = 6378  # Earth's radius (km)
wE = np.array([0, 0, 7.2921159e-5])  # Earth's angular velocity (rad/s)

# Satellite data
CD = 2.2
m = 100
A = np.pi / 4 * (1 ** 2)  # Frontal area (m^2)

# Initial orbital parameters (given)
rp = RE + 215  # Perigee radius (km)
ra = RE + 939  # Apogee radius (km)
RA = 339.94 * deg  # Right ascension of the node (radians)
i = 65.1 * deg  # Inclination (radians)
w = 58 * deg  # Argument of perigee (radians)
TA = 332 * deg  # True anomaly (radians)

# Initial orbital parameters (inferred)
e = (ra - rp) / (ra + rp)  # Eccentricity
a = (rp + ra) / 2  # Semimajor axis (km)
h = np.sqrt(mu * a * (1 - e ** 2))  # Angular momentum (km^2/s)
T = 2 * np.pi / np.sqrt(mu) * a ** 1.5  # Period (s)

# Store initial orbital elements
coe0 = [h, e, RA, i, w, TA]

# Get initial state vector
R0, V0 = sv_from_coe(coe0, mu)
r0 = np.linalg.norm(R0)  # Magnitude of R0
v0 = np.linalg.norm(V0)  # Magnitude of V0

# ODE45 to integrate the equations of motion
t0 = 0
tf = 120 * days
y0 = np.hstack((R0, V0))
nout = 40000
tspan = np.linspace(t0, tf, nout)

# Function to calculate the rates of change
def rates(t, f):
    R = f[:3]
    r = np.linalg.norm(R)
    alt = r - RE
    rho = atmosphere(alt)
    V = f[3:]
    Vrel = V - np.cross(wE, R)
    vrel = np.linalg.norm(Vrel)
    uv = Vrel / vrel
    ap = -CD * A / m * rho * (1000 * vrel) ** 2 / 2 * uv / 1000
    a0 = -mu * R / r ** 3
    a = a0 + ap
    return np.hstack((V, a))

# Function to specify the termination event
def terminate(t, y):
    R = y[:3]
    r = np.linalg.norm(R)
    alt = r - RE
    return alt - 100

terminate.terminal = True
terminate.direction = -1

# Solve the ODE
sol = solve_ivp(rates, [t0, tf], y0, method='RK45', t_eval=tspan, rtol=1e-8, atol=1e-8, events=terminate)

# Extract results
t = sol.t
y = sol.y.T
altitude = np.sqrt(np.sum(y[:, :3] ** 2, axis=1)) - RE

# Function to find extrema
def extrema(data):
    from scipy.signal import argrelextrema
    max_idx = argrelextrema(data, np.greater)[0]
    min_idx = argrelextrema(data, np.less)[0]
    return max_idx, min_idx

max_idx, min_idx = extrema(altitude)
max_altitude = altitude[max_idx]
min_altitude = altitude[min_idx]
maxima = np.vstack((t[max_idx], max_altitude)).T
minima = np.vstack((t[min_idx], min_altitude)).T
apogee = maxima[np.argsort(maxima[:, 0])]
perigee = minima[np.argsort(minima[:, 0])]

# Plot results
plt.figure()
plt.plot(apogee[:, 0] / days, apogee[:, 1], 'b', linewidth=2, label='Apogee')
plt.plot(perigee[:, 0] / days, perigee[:, 1], 'r', linewidth=2, label='Perigee')
plt.grid(True)
plt.xlabel('Time (days)')
plt.ylabel('Altitude (km)')
plt.ylim([0, 1000])
plt.legend()
plt.show()
