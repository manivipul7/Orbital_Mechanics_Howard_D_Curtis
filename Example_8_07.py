import numpy as np
from month_planet_names import month_planet_names
from planet_elements_and_sv import planet_elements_and_sv

# Global constants
mu = 1.327124e11  # gravitational parameter of the sun (km^3/s^2)
deg = np.pi / 180  # conversion factor between degrees and radians

# Input data
planet_id = 3
year = 2003
month = 8
day = 27
hour = 12
minute = 0
second = 0

# Algorithm 8.1: Calculate orbital elements and state vector
coe, r, v, jd = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second)

# Convert planet_id and month numbers into names for output (assuming month_planet_names exists)
month_name, planet_name = month_planet_names(month, planet_id)

# Print output as in MATLAB script
print('–––––––––––––––––––––––––––––––––––––––––––––––––––––')
print('\n Example 8.7')
print('\n\n Input data:\n')
print(f' Planet: {planet_name}')
print(f' Year : {year}')
print(f' Month : {month_name}')
print(f' Day : {day}')
print(f' Hour : {hour}')
print(f' Minute: {minute}')
print(f' Second: {second}')
print('\n\n Julian day: {:11.3f}'.format(jd))
print('\n\n Orbital elements:')
print('\n')
print(' Angular momentum (km^2/s) = {:g}'.format(coe[0]))
print(' Eccentricity = {:g}'.format(coe[1]))
print(' Right ascension of the ascending node (deg) = {:g}'.format(coe[2]))
print(' Inclination to the ecliptic (deg) = {:g}'.format(coe[3]))
print(' Argument of perihelion (deg) = {:g}'.format(coe[4]))
print(' True anomaly (deg) = {:g}'.format(coe[5]))
print(' Semimajor axis (km) = {:g}'.format(coe[6]))
print(' Longitude of perihelion (deg) = {:g}'.format(coe[7]))
print(' Mean longitude (deg) = {:g}'.format(coe[8]))
print(' Mean anomaly (deg) = {:g}'.format(coe[9]))
print(' Eccentric anomaly (deg) = {:g}'.format(coe[10]))
print('\n')
print(' State vector:')
print('\n')
print(' Position vector (km) = [{:g} {:g} {:g}]'.format(r[0], r[1], r[2]))
print(' Magnitude = {:g}'.format(np.linalg.norm(r)))
print('\n')
print(' Velocity (km/s)')
print(' = [{:g} {:g} {:g}]'.format(v[0], v[1], v[2]))
print(' Magnitude = {:g}'.format(np.linalg.norm(v)))
print('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')


