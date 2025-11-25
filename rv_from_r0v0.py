import numpy as np
from kepler_U import kepler_U
from f_and_g import f_and_g
from fDot_and_gDot import fDot_and_gDot

def rv_from_r0v0(R0, V0, t):
    """
    Computes the state vector (R, V) from the initial state vector (R0, V0) and the elapsed time.
    
    Parameters:
    R0 -- initial position vector (km)
    V0 -- initial velocity vector (km/s)
    t -- elapsed time (s)
    
    Returns:
    R -- final position vector (km)
    V -- final velocity vector (km/s)
    """
    global mu
    mu = 398600.0  # Gravitational parameter (km^3/s^2)
    
    # Magnitudes of R0 and V0
    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)
    
    # Initial radial velocity
    vr0 = np.dot(R0, V0) / r0
    
    # Reciprocal of the semimajor axis (from the energy equation)
    alpha = 2 / r0 - v0**2 / mu
    
    # Compute the universal anomaly
    x = kepler_U(t, r0, vr0, alpha, mu)
    
    # Compute the f and g functions
    f, g = f_and_g(x, t, r0, alpha, mu)
    
    # Compute the final position vector
    R = f * R0 + g * V0
    
    # Compute the magnitude of R
    r = np.linalg.norm(R)
    
    # Compute the derivatives of f and g
    fdot, gdot = fDot_and_gDot(x, r, r0, alpha, mu)
    
    # Compute the final velocity
    V = fdot * R0 + gdot * V0
    
    return R, V