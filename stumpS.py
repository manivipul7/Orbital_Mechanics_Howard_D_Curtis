import numpy as np

def stumpS(z):
    """
    Evaluates the Stumpff function S(z) according to Equation 3.52.
    
    Parameters:
    z : float
        Input argument.
        
    Returns:
    s : float
        Value of S(z).
    """
    
    if z > 0:
        s = (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z))**3
    elif z < 0:
        s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z))**3
    else:
        s = 1/6
    
    return s


