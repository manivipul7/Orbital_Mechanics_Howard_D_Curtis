import numpy as np
from sv_from_coe import sv_from_coe
from rva_relative import rva_relative

# Global constants
mu = 398600  # Gravitational parameter (km^3/s^2)
deg = np.pi / 180  # Degree to radian conversion

# Orbital parameters for Spacecraft A and B
h_A = 52059
e_A = 0.025724
i_A = 60 * deg
RAAN_A = 40 * deg
omega_A = 30 * deg
theta_A = 40 * deg

h_B = 52362
e_B = 0.0072696
i_B = 50 * deg
RAAN_B = 40 * deg
omega_B = 120 * deg
theta_B = 40 * deg

# Calculate state vectors for Spacecraft A and B
coe_A = [h_A, e_A, RAAN_A, i_A, omega_A, theta_A]
coe_B = [h_B, e_B, RAAN_B, i_B, omega_B, theta_B]

rA, vA = sv_from_coe(coe_A, mu)
rB, vB = sv_from_coe(coe_B, mu)

# Calculate relative position, velocity, and acceleration in LVLH frame
r, v, a = rva_relative(rA, vA, rB, vB, mu)

# Output results
print('\n\n––––––––––––––––––––––––––––––––––––––--–––––––––––––––––\n\n')
print('Orbital parameters of spacecraft A:')
print(f' angular momentum = {h_A} (km^2/s)')
print(f' eccentricity = {e_A}')
print(f' inclination = {i_A/deg} (deg)')
print(f' RAAN = {RAAN_A/deg} (deg)')
print(f' argument of perigee = {omega_A/deg} (deg)')
print(f' true anomaly = {theta_A/deg} (deg)\n')
print('State vector of spacecraft A:')
print(f' r = [{rA[0]}, {rA[1]}, {rA[2]}]')
print(f' (magnitude = {np.linalg.norm(rA)}) km')
print(f' v = [{vA[0]}, {vA[1]}, {vA[2]}]')
print(f' (magnitude = {np.linalg.norm(vA)}) km/s\n')

print('Orbital parameters of spacecraft B:')
print(f' angular momentum = {h_B} (km^2/s)')
print(f' eccentricity = {e_B}')
print(f' inclination = {i_B/deg} (deg)')
print(f' RAAN = {RAAN_B/deg} (deg)')
print(f' argument of perigee = {omega_B/deg} (deg)')
print(f' true anomaly = {theta_B/deg} (deg)\n')
print('State vector of spacecraft B:')
print(f' r = [{rB[0]}, {rB[1]}, {rB[2]}]')
print(f' (magnitude = {np.linalg.norm(rB)}) km')
print(f' v = [{vB[0]}, {vB[1]}, {vB[2]}]')
print(f' (magnitude = {np.linalg.norm(vB)}) km/s\n')

print('In the co-moving frame attached to A:')
print(f' Position of B relative to A = [{r[0]}, {r[1]}, {r[2]}]')
print(f' (magnitude = {np.linalg.norm(r)}) km\n')
print(f' Velocity of B relative to A = [{v[0]}, {v[1]}, {v[2]}]')
print(f' (magnitude = {np.linalg.norm(v)}) km/s\n')
print(f' Acceleration of B relative to A = [{a[0]}, {a[1]}, {a[2]}]')
print(f' (magnitude = {np.linalg.norm(a)}) km/s^2\n')

print('\n\n––––––––––––––––––––––––––––––––––––––--–––––––––––––––––\n\n')