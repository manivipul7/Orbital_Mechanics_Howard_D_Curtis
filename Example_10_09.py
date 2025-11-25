import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from solar_position import solar_position
from los import los
from sv_from_coe import sv_from_coe
from scipy.signal import savgol_filter

def rsmooth(y, window_length=101, polyorder=2):
    """Smooth data using a Savitzky-Golay filter."""
    return savgol_filter(y, window_length, polyorder)


# Constants
mu = 398600  # Gravitational parameter (km^3/s^2)
RE = 6378  # Earth's radius (km)
c = 2.998e8  # Speed of light (m/s)
S = 1367  # Solar constant (W/m^2)
Psr = S / c  # Solar pressure (Pa)
deg = np.pi / 180  # Degrees to radians conversion

# Satellite data
CR = 2  # Radiation pressure coefficient
m = 100  # Mass (kg)
As = 200  # Frontal area (m^2)

# Initial orbital parameters (given)
a0 = 10085.44  # Semimajor axis (km)
e0 = 0.025422  # Eccentricity
incl0 = np.deg2rad(88.3924)  # Inclination (radians)
RA0 = np.deg2rad(45.38124)  # Right ascension of the node (radians)
TA0 = np.deg2rad(343.4268)  # True anomaly (radians)
w0 = np.deg2rad(227.493)  # Argument of perigee (radians)

# Initial orbital parameters (inferred)
h0 = np.sqrt(mu * a0 * (1 - e0**2))  # Angular momentum (km^2/s)
T0 = 2 * np.pi / np.sqrt(mu) * a0**1.5  # Period (s)
rp0 = h0**2 / mu / (1 + e0)  # Perigee radius (km)
ra0 = h0**2 / mu / (1 - e0)  # Apogee radius (km)

# Store initial orbital elements in the vector coe0
coe0 = [h0, e0, RA0, incl0, w0, TA0]

# Time parameters
JD0 = 2438400.5  # Initial Julian date
t0 = 0  # Initial time (s)
tf = 3 * 365 * 24 * 3600  # Final time (s)
nout = 4000  # Number of solution points to output

# Integration time interval
tspan = np.linspace(t0, tf, nout)

# ODE options
reltol = 1e-8
abstol = 1e-8
initial_step = T0 / 1000

# Define the rates function
def rates(f, t):
    global JD
    JD = JD0 + t / (24 * 3600)  # Update the Julian Date at time t

    # Compute the apparent position vector of the sun
    lamda, eps, r_sun = solar_position(JD)
    lamda = np.deg2rad(lamda)  # Convert to radians
    eps = np.deg2rad(eps)  # Convert to radians

    # Extract the orbital elements at time t
    h, e, RA, incl, w, TA = f
    u = w + TA  # Argument of latitude

    # Compute the state vector at time t
    coe = [h, e, RA, incl, w, TA]
    R, V = sv_from_coe(coe, mu)

    # Calculate the magnitude of the radius vector
    r = np.linalg.norm(R)

    # Compute the shadow function and the solar radiation perturbation
    nu = los(R, r_sun)
    pSR = nu * (S / c) * CR * As / m / 1000

    # Calculate the trig functions in Equations 12.105
    sl = np.sin(lamda)
    cl = np.cos(lamda)
    se = np.sin(eps)
    ce = np.cos(eps)
    sW = np.sin(RA)
    si = np.sin(incl)
    su = np.sin(u)
    sT = np.sin(TA)
    cW = np.cos(RA)
    ci = np.cos(incl)
    cu = np.cos(u)
    cT = np.cos(TA)

    # Calculate the earth-sun unit vector components (Equations 12.105)
    ur = sl * ce * cW * ci * su + sl * ce * sW * cu - cl * sW * ci * su + cl * cW * cu + sl * se * si * su
    us = sl * ce * cW * ci * cu - sl * ce * sW * su - cl * sW * ci * cu - cl * cW * su + sl * se * si * cu
    uw = -sl * ce * cW * si + cl * sW * si + sl * se * ci

    # Calculate the time rates of the osculating elements from Equations 12.106
    hdot = -pSR * r * us
    edot = -pSR * (h / mu * sT * ur + 1 / mu / h * ((h**2 + mu * r) * cT + mu * e * r) * us)
    TAdot = h / r**2 - pSR / e / h * (h**2 / mu * cT * ur - (r + h**2 / mu) * sT * us)
    RAdot = -pSR * r / h / si * su * uw
    idot = -pSR * r / h * cu * uw
    wdot = -pSR * (-1 / e / h * (h**2 / mu * cT * ur - (r + h**2 / mu) * sT * us) - r * su / h / si * ci * uw)

    return [hdot, edot, RAdot, idot, wdot, TAdot]

# Initial conditions
y0 = coe0

# Integrate the equations using odeint
sol = odeint(rates, y0, tspan, rtol=reltol, atol=abstol, h0=initial_step)

# Extract the orbital elements' time histories from the solution vector y
h = sol[:, 0]
e = sol[:, 1]
RA = sol[:, 2]
incl = sol[:, 3]
w = sol[:, 4]
TA = sol[:, 5]
a = h**2 / mu / (1 - e**2)

# smooth the results
h = rsmooth(h)
e = rsmooth(e)
RA = rsmooth(RA)
incl = rsmooth(incl)
w = rsmooth(w)
TA = rsmooth(TA)
a = rsmooth(a)

# Plot the results
plt.figure(figsize=(10, 15))
plt.subplot(3, 2, 1)
plt.plot(tspan / (24 * 3600), h - h0)
plt.title('Angular Momentum (km^2/s)')
plt.xlabel('days')
plt.axis('tight')

plt.subplot(3, 2, 2)
plt.plot(tspan / (24 * 3600), e - e0)
plt.title('Eccentricity')
plt.xlabel('days')
plt.axis('tight')

plt.subplot(3, 2, 3)
plt.plot(tspan / (24 * 3600), a - a0)
plt.title('Semimajor axis (km)')
plt.xlabel('days')
plt.axis('tight')

plt.subplot(3, 2, 4)
plt.plot(tspan / (24 * 3600), (RA - RA0) / deg)
plt.title('Right Ascension (deg)')
plt.xlabel('days')
plt.axis('tight')

plt.subplot(3, 2, 5)
plt.plot(tspan / (24 * 3600), (incl - incl0) / deg)
plt.title('Inclination (deg)')
plt.xlabel('days')
plt.axis('tight')

plt.subplot(3, 2, 6)
plt.plot(tspan / (24 * 3600), (w - w0) / deg)
plt.title('Argument of Perigee (deg)')
plt.xlabel('days')
plt.axis('tight')

plt.show()
