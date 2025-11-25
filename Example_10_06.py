import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Conversion factors
hours = 3600
days = 24 * hours
deg = np.pi / 180

# Constants
mu = 398600  # Gravitational parameter (km^3/s^2)
RE = 6378  # Earth's radius (km)
J2 = 1082.63e-6  # Earth's J2

# Initial orbital parameters (given)
rp0 = RE + 300  # Perigee radius (km)
ra0 = RE + 3062  # Apogee radius (km)
RA0 = 45 * deg  # Right ascension of the node (radians)
i0 = 28 * deg  # Inclination (radians)
w0 = 30 * deg  # Argument of perigee (radians)
TA0 = 40 * deg  # True anomaly (radians)

# Initial orbital parameters (inferred)
e0 = (ra0 - rp0) / (ra0 + rp0)  # Eccentricity
h0 = np.sqrt(rp0 * mu * (1 + e0))  # Angular momentum (km^2/s)
a0 = (rp0 + ra0) / 2  # Semimajor axis (km)
T0 = 2 * np.pi / np.sqrt(mu) * a0 ** 1.5  # Period (s)

# Store initial orbital elements in the vector coe0
coe0 = [h0, e0, RA0, i0, w0, TA0]

# Time span for integration
t0 = 0
tf = 2 * days
nout = 5000  # Number of solution points to output for plotting purposes
tspan = np.linspace(t0, tf, nout)

# ODE options
options = {
    'rtol': 1e-8,
    'atol': 1e-8,
    'initial_step': T0 / 1000
}

# Define the rates function
def rates(t, f):
    h = f[0]
    e = f[1]
    RA = f[2]
    i = f[3]
    w = f[4]
    TA = f[5]
    
    r = h**2 / mu / (1 + e * np.cos(TA))  # The radius
    u = w + TA  # Argument of latitude
    
    # Orbital element rates at time t (Equations 10.89)
    hdot = -3/2 * J2 * mu * RE**2 / r**3 * np.sin(i)**2 * np.sin(2 * u)
    
    edot = 3/2 * J2 * mu * RE**2 / h / r**3 * (
        h**2 / mu / r * (np.sin(u) * np.sin(i)**2 * (3 * np.sin(TA) * np.sin(u) - 2 * np.cos(TA) * np.cos(u)) - np.sin(TA)) 
        - np.sin(i)**2 * np.sin(2 * u) * (e + np.cos(TA))
    )
    
    edot = 3/2 * J2 * mu * RE**2 / h / r**3 * (
        h**2 / mu / r * np.sin(TA) * (3 * np.sin(i)**2 * np.sin(u)**2 - 1)
        - np.sin(2 * u) * np.sin(i)**2 * ((2 + e * np.cos(TA)) * np.cos(TA) + e)
    )
    
    TAdot = h / r**2 + 3/2 * J2 * mu * RE**2 / e / h / r**3 * (
        h**2 / mu / r * np.cos(TA) * (3 * np.sin(i)**2 * np.sin(u)**2 - 1)
        + np.sin(2 * u) * np.sin(i)**2 * np.sin(TA) * (h**2 / mu / r + 1)
    )
    
    RAdot = -3 * J2 * mu * RE**2 / h / r**3 * np.sin(u)**2 * np.cos(i)
    
    idot = -3/4 * J2 * mu * RE**2 / h / r**3 * np.sin(2 * u) * np.sin(2 * i)
    
    wdot = 3/2 * J2 * mu * RE**2 / e / h / r**3 * (
        -h**2 / mu / r * np.cos(TA) * (3 * np.sin(i)**2 * np.sin(u)**2 - 1)
        - np.sin(2 * u) * np.sin(i)**2 * np.sin(TA) * (2 + e * np.cos(TA))
        + 2 * e * np.cos(i)**2 * np.sin(u)**2
    )
    
    # Pass these rates back to solve_ivp in the array dfdt
    dfdt = np.array([hdot, edot, RAdot, idot, wdot, TAdot])
    return dfdt

# Initial state
y0 = np.array(coe0)

# Integrate the Gauss variational equations
sol = solve_ivp(rates, [t0, tf], y0, t_eval=tspan, **options)

# Assign the time histories mnemonic variable names
t = sol.t
y = sol.y
h = y[0]
e = y[1]
RA = y[2]
i = y[3]
w = y[4]
TA = y[5]

# Plot the time histories of the osculating elements
plt.figure(1)
plt.subplot(5, 1, 1)
plt.plot(t / 3600, (RA - RA0) / deg)
plt.title('Right Ascension (degrees)')
plt.xlabel('hours')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.subplot(5, 1, 2)
plt.plot(t / 3600, (w - w0) / deg)
plt.title('Argument of Perigee (degrees)')
plt.xlabel('hours')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.subplot(5, 1, 3)
plt.plot(t / 3600, h - h0)
plt.title('Angular Momentum (km^2/s)')
plt.xlabel('hours')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.subplot(5, 1, 4)
plt.plot(t / 3600, e - e0)
plt.title('Eccentricity')
plt.xlabel('hours')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.subplot(5, 1, 5)
plt.plot(t / 3600, (i - i0) / deg)
plt.title('Inclination (degrees)')
plt.xlabel('hours')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.show()
