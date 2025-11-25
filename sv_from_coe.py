import numpy as np

def sv_from_coe(coe, mu):
    """
    This function computes the state vector (r, v) from the classical orbital elements (coe).

    Parameters:
    coe -- orbital elements [h, e, RA, incl, w, TA]
           where
           h   = angular momentum (km^2/s)
           e   = eccentricity
           RA  = right ascension of the ascending node (rad)
           incl= inclination of the orbit (rad)
           w   = argument of perigee (rad)
           TA  = true anomaly (rad)
    mu -- gravitational parameter (km^3/s^2)

    Returns:
    r -- position vector in the geocentric equatorial frame (km)
    v -- velocity vector in the geocentric equatorial frame (km/s)
    """
    h = coe[0]
    e = coe[1]
    RA = coe[2]
    incl = coe[3]
    w = coe[4]
    TA = coe[5]

    # Position and velocity vectors in the perifocal coordinate system
    rp = (h**2 / mu) * (1 / (1 + e * np.cos(TA))) * (np.array([np.cos(TA), np.sin(TA), 0]))
    vp = (mu / h) * (-np.sin(TA) * np.array([1, 0, 0]) + (e + np.cos(TA)) * np.array([0, 1, 0]))

    # Rotation matrices
    R3_W = np.array([
        [ np.cos(RA),  np.sin(RA), 0],
        [-np.sin(RA),  np.cos(RA), 0],
        [          0,           0, 1]
    ])
    
    R1_i = np.array([
        [1,           0,            0],
        [0,  np.cos(incl),  np.sin(incl)],
        [0, -np.sin(incl),  np.cos(incl)]
    ])

    R3_w = np.array([
        [ np.cos(w),  np.sin(w), 0],
        [-np.sin(w),  np.cos(w), 0],
        [         0,          0, 1]
    ])

    # Transformation matrix from perifocal to geocentric equatorial frame
    Q_pX = np.matmul(np.matmul(R3_w, R1_i), R3_W).T

    # Position and velocity vectors in the geocentric equatorial frame
    r = np.matmul(Q_pX, rp)
    v = np.matmul(Q_pX, vp)

    # Convert r and v into row vectors
    r = r.flatten()
    v = v.flatten()

    return r, v

