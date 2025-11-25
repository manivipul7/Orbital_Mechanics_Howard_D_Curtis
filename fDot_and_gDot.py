import numpy as np
from stumpC import stumpC
from stumpS import stumpS

def fDot_and_gDot(x, r, ro, a, mu):
    """
    Calculates the time derivatives of the Lagrange f and g coefficients.
    
    Parameters:
    x -- the universal anomaly after time t (km^0.5)
    r -- the radial position after time t (km)
    ro -- the radial position at time to (km)
    a -- reciprocal of the semimajor axis (1/km)
    
    Returns:
    fdot -- time derivative of the Lagrange f coefficient (1/s)
    gdot -- time derivative of the Lagrange g coefficient (dimensionless)
    """
    z = a * x**2
    
    # Equation 3.69c
    fdot = (np.sqrt(mu) / (r * ro)) * (z * stumpS(z) - 1) * x
    
    # Equation 3.69d
    gdot = 1 - (x**2 / r) * stumpC(z)
    
    return fdot, gdot