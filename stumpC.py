import numpy as np

def stumpC(z):
    """
    Evaluates the Stumpff function C(z) according to Equation 3.53.
    
    Parameters:
    z : float
        Input argument.
        
    Returns:
    c : float
        Value of C(z).
    """
    
    if z > 0:
        c = (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        c = 1/2
    
    return c


