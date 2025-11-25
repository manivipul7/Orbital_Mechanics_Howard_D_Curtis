import numpy as np

def f_and_g_ta(r0, v0, dt, mu):
    """
    This function calculates the Lagrange f and g coefficients from the change in true anomaly since time t0.

    Parameters:
    mu  -- gravitational parameter (km^3/s^2)
    dt  -- change in true anomaly (degrees)
    r0  -- position vector at time t0 (km)
    v0  -- velocity vector at time t0 (km/s)

    Returns:
    f -- the Lagrange f coefficient (dimensionless)
    g -- the Lagrange g coefficient (s)
    """
    h = np.linalg.norm(np.cross(r0, v0))
    vr0 = np.dot(v0, r0) / np.linalg.norm(r0)
    r0_mag = np.linalg.norm(r0)
    s = np.sin(np.radians(dt))
    c = np.cos(np.radians(dt))
    
    # Equation 2.152
    r = h**2 / mu / (1 + (h**2 / mu / r0_mag - 1) * c - h * vr0 * s / mu)
    
    # Equations 2.158a & b
    f = 1 - mu * r * (1 - c) / h**2
    g = r * r0_mag * s / h
    
    return f, g

# Example usage
if __name__ == '__main__':
    r0 = np.array([7000, 0, 0])  # Position vector at time t0 (km)
    v0 = np.array([0, 7.5, 1])   # Velocity vector at time t0 (km/s)
    dt = 30  # Change in true anomaly (degrees)
    mu = 398600  # Gravitational parameter (km^3/s^2)

    f, g = f_and_g_ta(r0, v0, dt, mu)
    print(f"f: {f}")
    print(f"g: {g}")
