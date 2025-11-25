import numpy as np
from numpy.linalg import norm
from scipy.optimize import root
from math import sin, cos, sqrt, pi, degrees
from gauss import gauss
from coe_from_sv import coe_from_sv

# Define constants
deg = pi / 180
mu = 398600  # Gravitational parameter (km^3/s^2)
Re = 6378  # Earth's radius (km)
f = 1 / 298.26  # Earth's flattening factor

# Function to convert degrees to radians
def deg2rad(deg):
    return deg * deg

# Function to convert radians to degrees
def rad2deg(rad):
    return rad / deg

# Example 5.11 data
H = 1
phi = 40 * deg
t = np.array([0, 118.104, 237.577])
ra = np.array([43.5365, 54.4196, 64.3178]) * deg
dec = np.array([-8.78334, -12.0739, -15.1054]) * deg
theta = np.array([44.5065, 45.000, 45.4992]) * deg

# Calculate site position vectors R and direction cosine vectors rho
fac1 = Re / sqrt(1 - (2 * f - f * f) * sin(phi) ** 2)
fac2 = (Re * (1 - f) ** 2 / sqrt(1 - (2 * f - f * f) * sin(phi) ** 2) + H) * sin(phi)
R = np.zeros((3, 3))
rho = np.zeros((3, 3))
for i in range(3):
    R[i, 0] = (fac1 + H) * cos(phi) * cos(theta[i])
    R[i, 1] = (fac1 + H) * cos(phi) * sin(theta[i])
    R[i, 2] = fac2
    rho[i, 0] = cos(dec[i]) * cos(ra[i])
    rho[i, 1] = cos(dec[i]) * sin(ra[i])
    rho[i, 2] = sin(dec[i])

# Run Gauss method to compute state vector
r, v, r_old, v_old = gauss(rho[0, :], rho[1, :], rho[2, :], R[0, :], R[1, :], R[2, :], t[0], t[1], t[2], mu)

# Calculate orbital elements from state vectors
coe_old = coe_from_sv(r_old, v_old, mu)
coe = coe_from_sv(r, v, mu)

# Print results
print("–––––––––––––––––––––––––––––––––––––––––––––––––––––")
print("Example 5.11: Orbit determination by the Gauss method")
print("\nRadius of earth (km):", Re)
print("Flattening factor:", f)
print("Gravitational parameter (km^3/s^2):", mu)
print("\nInput data:")
print("\nLatitude (deg):", rad2deg(phi))
print("Altitude above sea level (km):", H)
print("\nObservations:")
print("\nTime (s)\tAscension (deg)\tDeclination (deg)\tSidereal time (deg)")
for i in range(3):
    print(f"{t[i]:9.4f}\t{rad2deg(ra[i]):11.4f}\t{rad2deg(dec[i]):19.4f}\t{rad2deg(theta[i]):20.4f}")

print("\nSolution:")
print("\nWithout iterative improvement...")
print("\nr (km):", r_old)
print("v (km/s):", v_old)
print("\nAngular momentum (km^2/s):", coe_old[0])
print("Eccentricity:", coe_old[1])
print("RA of ascending node (deg):", rad2deg(coe_old[2]))
print("Inclination (deg):", rad2deg(coe_old[3]))
print("Argument of perigee (deg):", rad2deg(coe_old[4]))
print("True anomaly (deg):", rad2deg(coe_old[5]))
print("Semimajor axis (km):", coe_old[6])
print("Periapse radius (km):", coe_old[0] ** 2 / mu / (1 + coe_old[1]))
if coe_old[1] < 1:
    T = 2 * pi / sqrt(mu) * coe_old[6] ** 1.5
    print("Period:")
    print(f"Seconds: {T}")
    print(f"Minutes: {T / 60}")
    print(f"Hours: {T / 3600}")
    print(f"Days: {T / (24 * 3600)}")

print("\nWith iterative improvement...")
print("\nr (km):", r)
print("v (km/s):", v)
print("\nAngular momentum (km^2/s):", coe[0])
print("Eccentricity:", coe[1])
print("RA of ascending node (deg):", rad2deg(coe[2]))
print("Inclination (deg):", rad2deg(coe[3]))
print("Argument of perigee (deg):", rad2deg(coe[4]))
print("True anomaly (deg):", rad2deg(coe[5]))
print("Semimajor axis (km):", coe[6])
print("Periapase radius (km):", coe[0] ** 2 / mu / (1 + coe[1]))
if coe[1] < 1:
    T = 2 * pi / sqrt(mu) * coe[6] ** 1.5
    print("Period:")
    print(f"Seconds: {T}")
    print(f"Minutes: {T / 60}")
    print(f"Hours: {T / 3600}")
    print(f"Days: {T / (24 * 3600)}")

print("–––––––––––––––––––––––––––––––––––––––––––––––––––––")
