import numpy as np

def coe_from_sv(R, V, mu):
    """
    Computes the classical orbital elements (coe) from the state vector (R, V) using Algorithm 4.1.

    Parameters:
    R -- position vector in the geocentric equatorial frame (km)
    V -- velocity vector in the geocentric equatorial frame (km/s)
    mu -- gravitational parameter (km^3/s^2)

    Returns:
    coe -- vector of orbital elements [h, e, RA, incl, w, TA, a]
    """
    eps = 1.e-10
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)
    vr = np.dot(R, V) / r
    H = np.cross(R, V)
    h = np.linalg.norm(H)

    # Inclination (rad)
    incl = np.arccos(H[2] / h)

    # Node line (N)
    N = np.cross([0, 0, 1], H)
    n = np.linalg.norm(N)

    # Right ascension of the ascending node (RA) (rad)
    if n != 0:
        RA = np.arccos(N[0] / n)
        if N[1] < 0:
            RA = 2 * np.pi - RA
    else:
        RA = 0

    # Eccentricity vector (E) and eccentricity (e)
    E = (1 / mu) * ((v**2 - mu / r) * R - r * vr * V)
    e = np.linalg.norm(E)

    # Argument of perigee (w) (rad)
    if n != 0:
        if e > eps:
            w = np.arccos(np.dot(N, E) / (n * e))
            if E[2] < 0:
                w = 2 * np.pi - w
        else:
            w = 0
    else:
        w = 0

    # True anomaly (TA) (rad)
    if e > eps:
        TA = np.arccos(np.dot(E, R) / (e * r))
        if vr < 0:
            TA = 2 * np.pi - TA
    else:
        cp = np.cross(N, R)
        if cp[2] >= 0:
            TA = np.arccos(np.dot(N, R) / (n * r))
        else:
            TA = 2 * np.pi - np.arccos(np.dot(N, R) / (n * r))

    # Semi-major axis (a) (km)
    a = h**2 / mu / (1 - e**2)

    # Return orbital elements
    coe = [h, e, RA, incl, w, TA, a] 
    return coe


