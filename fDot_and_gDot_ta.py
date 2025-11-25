import numpy as np

def fDot_and_gDot_ta(r0, v0, dt, mu):
    """
    This function calculates the time derivatives of the Lagrange f and g coefficients from the change in true anomaly since time t0.

    Parameters:
    mu  -- gravitational parameter (km^3/s^2)
    dt  -- change in true anomaly (degrees)
    r0  -- position vector at time t0 (km)
    v0  -- velocity vector at time t0 (km/s)

    Returns:
    fdot -- time derivative of the Lagrange f coefficient (1/s)
    gdot -- time derivative of the Lagrange g coefficient (dimensionless)
    """
    h = np.linalg.norm(np.cross(r0, v0))
    vr0 = np.dot(v0, r0) / np.linalg.norm(r0)
    r0_mag = np.linalg.norm(r0)
    c = np.cos(np.radians(dt))
    s = np.sin(np.radians(dt))
    
    # Equations 2.158c & d
    fdot = mu / h * (vr0 / h * (1 - c) - s / r0_mag)
    gdot = 1 - mu * r0_mag / h**2 * (1 - c)
    
    return fdot, gdot

# Example usage
if __name__ == '__main__':
    r0 = np.array([7000, 0, 0])  # Position vector at time t0 (km)
    v0 = np.array([0, 7.5, 1])   # Velocity vector at time t0 (km/s)
    dt = 30  # Change in true anomaly (degrees)
    mu = 398600  # Gravitational parameter (km^3/s^2)

    fdot, gdot = fDot_and_gDot_ta(r0, v0, dt, mu)
    print(f"fdot: {fdot}")
    print(f"gdot: {gdot}")
