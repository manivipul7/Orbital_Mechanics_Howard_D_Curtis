import math

def kepler_H(e, M):
    """
    Solves Kepler's equation for the hyperbola e*sinh(F) - F = M using Newton's method.
    
    Parameters:
    e : float
        Eccentricity.
    M : float
        Hyperbolic mean anomaly (radians).
        
    Returns:
    F : float
        Hyperbolic eccentric anomaly (radians).
    """
    
    # Set an error tolerance
    error = 1.e-8
    
    # Starting value for F
    F = M
    
    # Iterate on Equation 3.45 until F is determined to within the error tolerance
    ratio = 1
    while abs(ratio) > error:
        ratio = (e * math.sinh(F) - F - M) / (e * math.cosh(F) - 1)
        F = F - ratio
    
    return F

