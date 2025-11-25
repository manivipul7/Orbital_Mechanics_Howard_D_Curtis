import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
from datetime import datetime, timedelta
from scipy.constants import pi, day
from simpsons_lunar_ephemeris import simpsons_lunar_ephemeris
from ra_and_dec_from_r import ra_and_dec_from_r

def juliandate(year, month, day, hour, minute, second):
    dt = datetime(year, month, day, hour, minute, second)
    jd = 367 * dt.year - int((7 * (dt.year + int((dt.month + 9) / 12))) / 4) + int((275 * dt.month) / 9) + dt.day + 1721013.5 + (dt.hour + dt.minute / 60 + dt.second / 3600) / 24
    return jd

def rates(t, y, jd0, days, ttt, mu_e, mu_m):
    jd = jd0 - (ttt - t) / days
    X, Y, Z, vX, vY, vZ = y
    r_ = np.array([X, Y, Z])
    r = np.linalg.norm(r_)
    rm_, _ = simpsons_lunar_ephemeris(jd)
    rms_ = rm_ - r_
    rms = np.linalg.norm(rms_)
    aearth_ = -mu_e * r_ / r**3
    amoon_ = mu_m * (rms_ / rms**3 - rm_ / np.linalg.norm(rm_)**3)
    a_ = aearth_ + amoon_
    dydt = [vX, vY, vZ, a_[0], a_[1], a_[2]]
    return dydt

def plotit_XYZ(X, Y, Z, Xm, Ym, Zm, imin):
    fig = plt.figure('Trajectories of Spacecraft (red) and Moon (green)', figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor((1, 1, 1))

    # Geocentric inertial coordinate axes
    L = 20 * Re
    ax.plot([0, L], [0, 0], [0, 0], color='k')
    ax.text(L, 0, 0, 'X', fontsize=12, fontstyle='italic', fontname='Palatino')
    ax.plot([0, 0], [0, L], [0, 0], color='k')
    ax.text(0, L, 0, 'Y', fontsize=12, fontstyle='italic', fontname='Palatino')
    ax.plot([0, 0], [0, 0], [0, L], color='k')
    ax.text(0, 0, L, 'Z', fontsize=12, fontstyle='italic', fontname='Palatino')

    # Earth
    u, v = np.mgrid[0:2*np.pi:128j, 0:np.pi:64j]
    xx = Re * np.cos(u) * np.sin(v)
    yy = Re * np.sin(u) * np.sin(v)
    zz = Re * np.cos(v)
    ax.plot_surface(xx, yy, zz, color='b', alpha=0.5)

    # Spacecraft at TLI
    ax.plot([X[0]], [Y[0]], [Z[0]], 'o', markeredgecolor='k', markerfacecolor='k', markersize=3)

    # Spacecraft at closest approach
    ax.plot([X[imin]], [Y[imin]], [Z[imin]], 'o', markeredgecolor='k', markerfacecolor='k', markersize=2)

    # Spacecraft at tf
    ax.plot([X[-1]], [Y[-1]], [Z[-1]], 'o', markeredgecolor='r', markerfacecolor='r', markersize=3)

    # Moon at TLI
    ax.text(Xm[0], Ym[0], Zm[0], 'Moon at TLI')
    xx = Rm * np.cos(u) * np.sin(v)
    yy = Rm * np.sin(u) * np.sin(v)
    zz = Rm * np.cos(v)
    ax.plot_surface(xx + Xm[0], yy + Ym[0], zz + Zm[0], color='g', alpha=0.99)

    # Moon at closest approach
    ax.plot_surface(xx + Xm[imin], yy + Ym[imin], zz + Zm[imin], color='g', alpha=0.99)

    # Moon at end of simulation
    ax.plot_surface(xx + Xm[-1], yy + Ym[-1], zz + Zm[-1], color='g', alpha=0.99)

    # Spacecraft trajectory
    ax.plot(X, Y, Z, 'r', linewidth=1.5)

    # Moon trajectory
    ax.plot(Xm, Ym, Zm, 'g', linewidth=0.5)

    ax.set_aspect('auto')
    ax.axis('off')
    ax.view_init(elev=20, azim=30)
    plt.show()
    

def plotit_xyz(x, y, z, xm, ym, zm, imin):
    # Placeholder function for plotting the trajectory in moon-fixed frame
    pass

# Constants
deg = pi / 180
days = 24 * 3600
Re = 6378
Rm = 1737
m_e = 5974.e21
m_m = 73.48e21
mu_e = 398600.4
mu_m = 4902.8
D = 384400
RS = D * (m_m / m_e)**(2 / 5)

# Data for Example 9.03
year = 2020
month = 5
day = 4
hour = 12
minute = 0
second = 0
t0 = 0
z0 = 320
alpha0 = 90
dec0 = 15
gamma0 = 40
fac = .9924  # Fraction of Vesc
ttt = 3 * days
tf = ttt + 2.667 * days

# State vector of moon at target date:
jd0 = juliandate(year, month, day, hour, minute, second)
rm0_, vm0_ = simpsons_lunar_ephemeris(jd0)
RA, Dec = ra_and_dec_from_r(rm0_)
distance = np.linalg.norm(rm0_)
hmoon_ = np.cross(rm0_, vm0_)
hmoon = np.linalg.norm(hmoon_)
inclmoon = np.arccos(hmoon_[2] / hmoon) * 180 / pi

# Initial position vector of probe:
I_ = np.array([1, 0, 0])
J_ = np.array([0, 1, 0])
K_ = np.cross(I_, J_)
r0 = Re + z0
r0_ = r0 * (np.cos(alpha0 * deg) * np.cos(dec0 * deg) * I_ +
            np.sin(alpha0 * deg) * np.cos(dec0 * deg) * J_ +
            np.sin(dec0 * deg) * K_)
vesc = np.sqrt(2 * mu_e / r0)
v0 = fac * vesc
w0_ = np.cross(r0_, rm0_) / np.linalg.norm(np.cross(r0_, rm0_))

# Initial velocity vector of probe:
ur_ = r0_ / np.linalg.norm(r0_)
uperp_ = np.cross(w0_, ur_) / np.linalg.norm(np.cross(w0_, ur_))
vr = v0 * np.sin(gamma0 * deg)
vperp = v0 * np.cos(gamma0 * deg)
v0_ = vr * ur_ + vperp * uperp_
uv0_ = v0_ / v0

# Initial state vector of the probe:
y0 = np.concatenate((r0_, v0_))

# Solving the ODE
solution = solve_ivp(rates, [t0, tf], y0, args=(jd0, days, ttt, mu_e, mu_m), rtol=1e-10, atol=1e-10)
t = solution.t
y = solution.y.T

# Extracting positions and velocities
X, Y, Z, vX, vY, vZ = y[:, 0], y[:, 1], y[:, 2], y[:, 3], y[:, 4], y[:, 5]

# Moon and probe trajectories (initializing arrays)
xm, ym, zm = [], [], []
Xm, Ym, Zm = [], [], []
vXm, vYm, vZm = [], [], []

# Finding perilune
dist_min = 1e30
for i, ti in enumerate(t):
    r_ = np.array([X[i], Y[i], Z[i]])
    jd = jd0 - (ttt - ti) / days
    rm_, vm_ = simpsons_lunar_ephemeris(jd)
    Xm.append(rm_[0])
    Ym.append(rm_[1])
    Zm.append(rm_[2])
    vXm.append(vm_[0])
    vYm.append(vm_[1])
    vZm.append(vm_[2])

    x_ = rm_
    z_ = np.cross(x_, vm_)
    y_ = np.cross(z_, x_)
    i_ = x_ / np.linalg.norm(x_)
    j_ = y_ / np.linalg.norm(y_)
    k_ = z_ / np.linalg.norm(z_)
    Q = np.vstack([i_, j_, k_])
    rx_ = Q @ r_
    xm.append(rx_[0])
    ym.append(rx_[1])
    zm.append(rx_[2])
    rmx_ = Q @ rm_
    dist_ = r_ - rm_
    dist = np.linalg.norm(dist_)
    if dist < dist_min:
        imin = i
        dist_min_ = dist_
        dist_min = dist

rmTLI_ = np.array([Xm[0], Ym[0], Zm[0]])
RATLI, DecTLI = ra_and_dec_from_r(rmTLI_)
v_atdmin_ = np.array([vX[imin], vY[imin], vZ[imin]])
rm_perilune_ = np.array([Xm[imin], Ym[imin], Zm[imin]])
vm_perilune_ = np.array([vXm[imin], vYm[imin], vZm[imin]])
RA_at_perilune, Dec_at_perilune = ra_and_dec_from_r(rm_perilune_)
target_error = np.linalg.norm(rm_perilune_ - rm0_)
rel_speed = np.linalg.norm(v_atdmin_ - vm_perilune_)
rend_ = np.array([X[-1], Y[-1], Z[-1]])
alt_end = np.linalg.norm(rend_) - Re
ra_end, dec_end = ra_and_dec_from_r(rend_)

# Output results
print(f"\n\nExample 9.3 4e\n")
print(f"Date and time of arrival at moon: {month}/{day}/{year} {hour}:{minute}:{second}")
print(f"Moon's position:")
print(f" Distance = {distance:11g} km")
print(f" Right Ascension = {RA:11g} deg")
print(f" Declination = {Dec:11g} deg")
print(f" Moon's orbital inclination = {inclmoon:11g} deg\n")

print(f"The probe at earth departure (t = {t0} sec):")
print(f" Altitude = {z0:11g} km")
print(f" Right ascension = {alpha0:11g} deg")
print(f" Declination = {dec0:11g} deg")
print(f" Flight path angle = {gamma0:11g} deg")
print(f" Speed = {v0:11g} km/s")
print(f" Escape speed = {vesc:11g} km/s")
print(f" v/vesc = {v0/vesc}")
print(f" Inclination of translunar orbit = {np.degrees(np.arccos(w0_[2])):11g} deg\n")

print(f"The moon when the probe is at TLI:")
print(f" Distance = {np.linalg.norm(rmTLI_):11g} km")
print(f" Right Ascension = {RATLI:11g} deg")
print(f" Declination = {DecTLI:11g} deg\n")

print(f"The moon when the probe is at perilune: ")
print(f" Distance = {np.linalg.norm(rm_perilune_):11g} km")
print(f" Speed = {np.linalg.norm(vm_perilune_):11g} km/s")
print(f" Right Ascension = {RA_at_perilune:11g} deg")
print(f" Declination = {Dec_at_perilune:11g} deg")
print(f" Target error = {target_error:11g} km\n\n")

print(f"The probe at perilune:")
print(f" Altitude = {dist_min - Rm:11g} km")
print(f" Speed = {np.linalg.norm(v_atdmin_):11g} km/s")
print(f" Relative speed = {rel_speed:11g} km/s\n")
print(f"Inclination of osculating plane = {incl[imin]:11g} deg")
print(f" Time from TLI to perilune = {abs(t[imin])/3600:11g} hours ({abs(t[imin])/3600/24} days)")

print(f"\n\nTotal time of flight = {t[-1]/days:11g} days")
print(f"Time to target point = {ttt/days:11g} days")
print(f"Final earth altitude = {alt_end:11g} km")
print(f"Final right ascension = {ra_end:11g} deg")
print(f"Final declination = {dec_end:11g} deg\n")

plt.rcParams['font.family'] = 'Arial'

# Function to plot the trajectory in the inertial frame
# def plotit_XYZ(X, Y, Z, Xm, Ym, Zm, imin):
#     fig = plt.figure('Trajectories of Spacecraft (red) and Moon (green)')
#     ax = fig.add_subplot(111, projection='3d')
#     ax.set_facecolor((1, 1, 1))
# 
#     u, v = np.mgrid[0:2*np.pi:128j, 0:np.pi:128j]
#     xx = np.cos(u) * np.sin(v)
#     yy = np.sin(u) * np.sin(v)
#     zz = np.cos(v)
# 
#     # Geocentric inertial coordinate axes:
#     L = 20 * Re
#     ax.plot([0, L], [0, 0], [0, 0], color='k')
#     ax.text(L, 0, 0, 'X', fontsize=12, fontstyle='italic', fontname='Arial')
#     ax.plot([0, 0], [0, L], [0, 0], color='k')
#     ax.text(0, L, 0, 'Y', fontsize=12, fontstyle='italic', fontname='Arial')
#     ax.plot([0, 0], [0, 0], [0, L], color='k')
#     ax.text(0, 0, L, 'Z', fontsize=12, fontstyle='italic', fontname='Arial')
# 
#     # Earth
#     ax.plot_surface(Re * xx, Re * yy, Re * zz, color='b', alpha=0.5, shade=True)
# 
#     # Spacecraft at TLI
#     ax.plot([X[0]], [Y[0]], [Z[0]], 'ko', markersize=3)
#     # Spacecraft at closest approach
#     ax.plot([X[imin]], [Y[imin]], [Z[imin]], 'ko', markersize=2)
#     # Spacecraft at tf
#     ax.plot([X[-1]], [Y[-1]], [Z[-1]], 'ro', markersize=3)
# 
#     # Moon at TLI
#     ax.text(Xm[0], Ym[0], Zm[0], 'Moon at TLI')
#     ax.plot_surface(Rm * xx + Xm[0], Rm * yy + Ym[0], Rm * zz + Zm[0], color='g', alpha=0.99, shade=True)
# 
#     # Moon at closest approach
#     ax.plot_surface(Rm * xx + Xm[imin], Rm * yy + Ym[imin], Rm * zz + Zm[imin], color='g', alpha=0.99, shade=True)
# 
#     # Moon at end of simulation
#     ax.plot_surface(Rm * xx + Xm[-1], Rm * yy + Ym[-1], Rm * zz + Zm[-1], color='g', alpha=0.99, shade=True)
# 
#     # Spacecraft trajectory
#     ax.plot(X, Y, Z, 'r', linewidth=1.5)
# 
#     # Moon trajectory
#     ax.plot(Xm, Ym, Zm, 'g', linewidth=0.5)
# 
#     ax.set_aspect('auto')
#     ax.set_xlabel('X (km)')
#     ax.set_ylabel('Y (km)')
#     ax.set_zlabel('Z (km)')
#     plt.show()

# Function to plot the trajectory in the Moon-fixed rotating frame
def plotit_xyz(x, y, z, xm, ym, zm, imin):
    fig = plt.figure('Spacecraft trajectory in Moon-fixed rotating frame')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor((1, 1, 1))

    u, v = np.mgrid[0:2*np.pi:128j, 0:np.pi:128j]
    xx = np.cos(u) * np.sin(v)
    yy = np.sin(u) * np.sin(v)
    zz = np.cos(v)

    # Spacecraft trajectory
    ax.plot(x, y, z, 'r', linewidth=2.0)

    # Moon trajectory
    ax.plot(xm, ym, zm, 'g', linewidth=0.5)

    # Earth
    ax.plot_surface(Re * xx, Re * yy, Re * zz, color='b', alpha=0.5, shade=True)

    # Geocentric moon-fixed coordinate axes
    L1 = 63 * Re; L2 = 20 * Re; L3 = 29 * Re
    ax.plot([0, L1], [0, 0], [0, 0], 'k')
    ax.text(L1, 0, 0, 'x', fontsize=12, fontstyle='italic', fontname='Arial')
    ax.plot([0, 0], [0, L2], [0, 0], 'k')
    ax.text(0, L2, 0, 'y', fontsize=12, fontstyle='italic', fontname='Arial')
    ax.plot([0, 0], [0, 0], [0, L3], 'k')
    ax.text(0, 0, L3, 'z', fontsize=12, fontstyle='italic', fontname='Arial')

    # Spacecraft at TLI
    ax.plot([x[0]], [y[0]], [z[0]], 'ko', markersize=3)
    # Spacecraft at closest approach
    ax.plot([x[imin]], [y[imin]], [z[imin]], 'ko', markersize=2)
    # Spacecraft at tf
    ax.plot([x[-1]], [y[-1]], [z[-1]], 'ro', markersize=3)

    # Moon at TLI
    ax.text(xm[0], ym[0], zm[0], 'Moon at TLI')
    ax.plot_surface(Rm * xx + xm[0], Rm * yy + ym[0], Rm * zz + zm[0], color='g', alpha=0.99, shade=True)

    # Moon at spacecraft closest approach
    ax.plot_surface(Rm * xx + xm[imin], Rm * yy + ym[imin], Rm * zz + zm[imin], color='g', alpha=0.99, shade=True)

    # Moon at end of simulation
    ax.plot_surface(Rm * xx + xm[-1], Rm * yy + ym[-1], Rm * zz + zm[-1], color='g', alpha=0.99, shade=True)

    ax.set_aspect('auto')
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_zlabel('z (km)')
    plt.show()

# Function to plot inclination vs distance from Moon
def plot_inclination_vs_distance(rms, incl):
    plt.figure()
    plt.plot(rms / RS, incl)
    plt.axhline(y=90, color='r', linestyle='-')
    plt.title('Osculating Plane Inclination vs Distance from Moon')
    plt.xlabel(r'$r_{ms}/R_s$')
    plt.ylabel('Inclination (deg)')
    plt.grid(True)
    plt.grid(which='minor')
    plt.show()

# Example usage
# Replace with actual data
X = np.array([0, 1, 2])
Y = np.array([0, 1, 2])
Z = np.array([0, 1, 2])
Xm = np.array([0, 1, 2])
Ym = np.array([0, 1, 2])
Zm = np.array([0, 1, 2])
imin = 1

x = np.array([0, 1, 2])
y = np.array([0, 1, 2])
z = np.array([0, 1, 2])
xm = np.array([0, 1, 2])
ym = np.array([0, 1, 2])
zm = np.array([0, 1, 2])

rms = np.array([0, 1, 2])
incl = np.array([0, 45, 90])

plotit_XYZ(X, Y, Z, Xm, Ym, Zm, imin)
plot_inclination_vs_distance(rms, incl)
plotit_xyz(x, y, z, xm, ym, zm, imin)



