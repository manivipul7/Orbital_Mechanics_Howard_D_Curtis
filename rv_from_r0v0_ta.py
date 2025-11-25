# rv_from_r0v0_ta.py
import numpy as np
from f_and_g_ta import f_and_g_ta
from fDot_and_gDot_ta import fDot_and_gDot_ta

def rv_from_r0v0_ta(r0, v0, dt, mu):
    """
    This function computes the state vector (r, v) from the
    initial state vector (r0, v0) and the change in true anomaly.

    Parameters:
    mu  -- gravitational parameter (km^3/s^2)
    r0  -- initial position vector (km)
    v0  -- initial velocity vector (km/s)
    dt  -- change in true anomaly (degrees)

    Returns:
    r -- final position vector (km)
    v -- final velocity vector (km/s)
    """
    # Compute the f and g functions and their derivatives
    f, g = f_and_g_ta(r0, v0, dt, mu)
    fdot, gdot = fDot_and_gDot_ta(r0, v0, dt, mu)

    # Compute the final position and velocity vectors
    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0
    
    return r, v

# Example usage
if __name__ == '__main__':
    r0 = np.array([7000, 0, 0])  # Initial position vector (km)
    v0 = np.array([0, 7.5, 1])   # Initial velocity vector (km/s)
    dt = 30  # Change in true anomaly (degrees)
    mu = 398600  # Gravitational parameter (km^3/s^2)

    r, v = rv_from_r0v0_ta(r0, v0, dt, mu)
    print(f"Final position vector: {r}")
    print(f"Final velocity vector: {v}")
