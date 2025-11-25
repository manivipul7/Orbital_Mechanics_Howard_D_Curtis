import numpy as np

def ra_and_dec_from_r(r):
    """
    Calculates the right ascension and the declination from the geocentric equatorial position vector.

    Parameters:
    r -- position vector (numpy array or list of length 3)

    Returns:
    ra -- right ascension (degrees)
    dec -- declination (degrees)
    """
    r_norm = np.linalg.norm(r)
    l = r[0] / r_norm
    m = r[1] / r_norm
    n = r[2] / r_norm

    dec = np.degrees(np.arcsin(n))

    if m > 0:
        ra = np.degrees(np.arccos(l / np.cos(np.radians(dec))))
    else:
        ra = 360 - np.degrees(np.arccos(l / np.cos(np.radians(dec))))

    return ra, dec