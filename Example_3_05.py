import math
from kepler_H import kepler_H

# Example data from Example 3.5
e = 2.7696
M = 40.69

# Compute hyperbolic eccentric anomaly F using kepler_H function
F = kepler_H(e, M)

# Print results to the console
print(f"For eccentricity e = {e} and hyperbolic mean anomaly M = {M}:")
print(f"Hyperbolic eccentric anomaly F = {F}")