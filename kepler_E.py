import math

def kepler_E(e, M):
    """
    Solves Kepler's equation E - e*sin(E) = M using Newton's method.
    
    Parameters:
    e : float
        Eccentricity.
    M : float
        Mean anomaly (radians).
        
    Returns:
    E : float
        Eccentric anomaly (radians).
    """
    
    # Set an error tolerance
    error = 1.e-8
    
    # Select a starting value for E
    if M < math.pi:
        E = M + e / 2
    else:
        E = M - e / 2
    
    # Iterate on Equation 3.17 until E is determined to within the error tolerance
    ratio = 1
    while abs(ratio) > error:
        ratio = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
        E = E - ratio
    
    return E
