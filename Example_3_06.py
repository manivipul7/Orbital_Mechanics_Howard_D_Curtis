from kepler_U import kepler_U
import math
# Define global variable mu (gravitational parameter), assumed to be defined elsewhere
global mu
# Example data for Example_3_06
dt = 3600.0  # Time since x = 0 (s)
ro = 10000.0  # Radial position when x = 0 (km)
vro = 3.0752  # Radial velocity when x = 0 (km/s)
a = -19655.0  # Semimajor axis (km)
# Call kepler_U function to solve for the universal anomaly x
mu = 398600.0  # Gravitational parameter (km^3/s^2)
x = kepler_U(dt, ro, vro, 1.0 / a, mu)


# Output results
print("\n Example 3.6")
print(f" Initial radial coordinate (km) = {ro}")
print(f" Initial radial velocity (km/s) = {vro}")
print(f" Elapsed time (seconds) = {dt}")
print(f" Semimajor axis (km) = {a}")
print(f" Universal anomaly (km^0.5) = {x}")
print("\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n")