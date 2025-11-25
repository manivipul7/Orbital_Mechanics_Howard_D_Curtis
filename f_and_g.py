import numpy as np
from stumpC import stumpC
from stumpS import stumpS

def f_and_g(x, t, ro, a, mu):
    """
    Calculates the Lagrange f and g coefficients.
    
    Parameters:
    x -- the universal anomaly after time t (km^0.5)
    t -- the time elapsed since ro (s)
    ro -- the radial position at time to (km)
    a -- reciprocal of the semimajor axis (1/km)
    
    Returns:
    f -- the Lagrange f coefficient (dimensionless)
    g -- the Lagrange g coefficient (s)
    """
    z = a * x**2
    
    # Equation 3.69a
    f = 1 - (x**2 / ro) * stumpC(z)
    
    # Equation 3.69b
    g = t - (1 / np.sqrt(mu)) * (x**3) * stumpS(z)
    
    return f, g