import numpy as np
from coe_from_sv import coe_from_sv

# Constants and initial data
deg = np.pi / 180
mu = 398600

# Data declaration for Example 4.3
r = np.array([-6045, -3490, 2500])
v = np.array([-3.457, 6.618, 2.533])

# Compute orbital elements
coe = coe_from_sv(r, v, mu)

# Unpack orbital elements
h, e, RA, incl, w, TA, a = coe

# Print results
print('–––––––––––––––––––––––––––––––––––––––––––––––––––––')
print('\n Example 4.3\n')
print(f'\n Gravitational parameter (km^3/s^2) = {mu}')
print('\n State vector:')
print(f'\n r (km) = [{r[0]} {r[1]} {r[2]}]')
print(f'\n v (km/s) = [{v[0]} {v[1]} {v[2]}]')
print('\n Angular momentum (km^2/s) =', h)
print('\n Eccentricity =', e)
print('\n Right ascension (deg) =', np.degrees(RA))
print('\n Inclination (deg) =', np.degrees(incl))
print('\n Argument of perigee (deg) =', np.degrees(w))
print('\n True anomaly (deg) =', np.degrees(TA))
print('\n Semimajor axis (km) =', a)
print('\n')

# Calculate and print the period if the orbit is an ellipse
if e < 1:
    T = 2 * np.pi / np.sqrt(mu) * a**1.5
    print('\n Period:')
    print(f'\n    Seconds = {T}')
    print(f'\n    Minutes = {T / 60}')
    print(f'\n    Hours = {T / 3600}')
    print(f'\n    Days = {T / 24 / 3600}')
print('–––––––––––––––––––––––––––––––––––––––––––––––––––––')