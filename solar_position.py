import numpy as np

def solar_position(jd):
    """
    This function calculates the geocentric equatorial position vector
    of the sun, given the Julian date.

    Returns:
        lamda (float): Apparent ecliptic longitude (degrees)
        eps (float): Obliquity of the ecliptic (degrees)
        r_S (numpy array): Geocentric position vector (km)
    """
    # Astronomical unit (km)
    AU = 149597870.691
    
    # Julian days since J2000
    n = jd - 2451545
    
    # Julian centuries since J2000
    cy = n / 36525
    
    # Mean anomaly (degrees)
    M = 357.528 + 0.9856003 * n
    M = M % 360
    
    # Mean longitude (degrees)
    L = 280.460 + 0.98564736 * n
    L = L % 360
    
    # Apparent ecliptic longitude (degrees)
    lamda = L + 1.915 * np.sin(np.radians(M)) + 0.020 * np.sin(np.radians(2 * M))
    lamda = lamda % 360
    
    # Obliquity of the ecliptic (degrees)
    eps = 23.439 - 0.0000004 * n
    
    # Unit vector from earth to sun
    u = np.array([
        np.cos(np.radians(lamda)), 
        np.sin(np.radians(lamda)) * np.cos(np.radians(eps)), 
        np.sin(np.radians(lamda)) * np.sin(np.radians(eps))
    ])
    
    # Distance from earth to sun (km)
    rS = (1.00014 - 0.01671 * np.cos(np.radians(M)) - 0.000140 * np.cos(np.radians(2 * M))) * AU
    
    # Geocentric position vector (km)
    r_S = rS * u
    
    return lamda, eps, r_S
