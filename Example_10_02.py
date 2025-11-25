import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from sv_from_coe import sv_from_coe
from rv_from_r0v0 import rv_from_r0v0
from coe_from_sv import coe_from_sv

# Global constants
mu = 398600  # Gravitational parameter (km^3/s^2)
RE = 6378  # Earth's radius (km)
J2 = 1082.63e-6  # J2 perturbation coefficient

# Conversion factors
hours = 3600
days = 24 * hours
deg = np.pi / 180

# Initial orbital parameters (given)
zp0 = 300  # Perigee altitude (km)
za0 = 3062  # Apogee altitude (km)
RA0 = 45 * deg  # Right ascension of the node (radians)
i0 = 28 * deg  # Inclination (radians)
w0 = 30 * deg  # Argument of perigee (radians)
TA0 = 40 * deg  # True anomaly (radians)

# Initial orbital parameters (inferred)
rp0 = RE + zp0  # Perigee radius (km)
ra0 = RE + za0  # Apogee radius (km)
e0 = (ra0 - rp0) / (ra0 + rp0)  # Eccentricity
a0 = (ra0 + rp0) / 2  # Semimajor axis (km)
h0 = np.sqrt(rp0 * mu * (1 + e0))  # Angular momentum (km^2/s)
T0 = 2 * np.pi / np.sqrt(mu) * a0 ** 1.5  # Period (s)

t0 = 0
tf = 2 * days

# Store the initial orbital elements
coe0 = [h0, e0, RA0, i0, w0, TA0]

# Get initial state vector
R0, V0 = sv_from_coe(coe0, mu)
r0 = np.linalg.norm(R0)  # Magnitude of R0
v0 = np.linalg.norm(V0)  # Magnitude of V0

del_t = T0 / 100  # Time step for Encke procedure

# ODE settings
def rates(t, f):
    del_r = f[:3]
    del_v = f[3:]
    
    # Compute the state vector on the osculating orbit at time t
    Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0)
    
    # Calculate the components of the state vector on the perturbed orbit
    Rpp = Rosc + del_r
    Vpp = Vosc + del_v
    rosc = np.linalg.norm(Rosc)
    rpp = np.linalg.norm(Rpp)
    
    # Compute the J2 perturbing acceleration
    xx, yy, zz = Rpp
    fac = 3/2 * J2 * (mu / rpp**2) * (RE / rpp)**2
    ap = -fac * np.array([(1 - 5*(zz/rpp)**2)*(xx/rpp), 
                          (1 - 5*(zz/rpp)**2)*(yy/rpp), 
                          (3 - 5*(zz/rpp)**2)*(zz/rpp)])
    
    # Compute the total perturbing acceleration
    F = 1 - (rosc/rpp)**3
    del_a = -mu / rosc**3 * (del_r - F * Rpp) + ap
    
    return np.hstack((del_v, del_a))

# Encke integration
t0 = 0
tf = 2 * days
t = t0
tsave = [t0]
y = np.hstack((R0, V0))
del_y0 = np.zeros(6)
t += del_t

# Time integration loop
while t <= tf + del_t / 2:
    sol = solve_ivp(rates, [t0, t], del_y0, method='RK45', max_step=del_t)
    z = sol.y[:, -1]
    
    Rosc, Vosc = rv_from_r0v0(R0, V0, t - t0)
    R0 = Rosc + z[:3]
    V0 = Vosc + z[3:]
    t0 = t
    
    tsave.append(t)
    y = np.vstack((y, np.hstack((R0, V0))))
    
    t += del_t
    del_y0 = np.zeros(6)

t = np.array(tsave)

# Extracting orbital elements from the state vector at each solution time
n_times = len(t)
r = np.zeros(n_times)
v = np.zeros(n_times)
h = np.zeros(n_times)
e = np.zeros(n_times)
RA = np.zeros(n_times)
i = np.zeros(n_times)
w = np.zeros(n_times)
TA = np.zeros(n_times)
a = np.zeros(n_times)

for j in range(n_times):
    R = y[j, :3]
    V = y[j, 3:]
    r[j] = np.linalg.norm(R)
    v[j] = np.linalg.norm(V)
    coe = coe_from_sv(R, V, mu)
    h[j], e[j], RA[j], i[j], w[j], TA[j], a[j]= coe

# Plotting selected osculating elements
plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(t / 3600, (RA - RA0) / deg)
plt.title('Variation of Right Ascension')
plt.xlabel('hours')
plt.ylabel(r'$\Delta\Omega$ (deg)')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.subplot(2, 1, 2)
plt.plot(t / 3600, (w - w0) / deg)
plt.title('Variation of Argument of Perigee')
plt.xlabel('hours')
plt.ylabel(r'$\Delta\omega$ (deg)')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.figure(2)
plt.subplot(3, 1, 1)
plt.plot(t / 3600, h - h0)
plt.title('Variation of Angular Momentum')
plt.xlabel('hours')
plt.ylabel(r'$\Delta h$ (km$^2$/s)')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.subplot(3, 1, 2)
plt.plot(t / 3600, e - e0)
plt.title('Variation of Eccentricity')
plt.xlabel('hours')
plt.ylabel(r'$\Delta e$')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.subplot(3, 1, 3)
plt.plot(t / 3600, (i - i0) / deg)
plt.title('Variation of Inclination')
plt.xlabel('hours')
plt.ylabel(r'$\Delta i$ (deg)')
plt.grid(True)
plt.grid(which='minor')
plt.axis('tight')

plt.show()
