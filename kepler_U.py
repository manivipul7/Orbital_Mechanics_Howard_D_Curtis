import math
from stumpC import stumpC
from stumpS import stumpS

def kepler_U(dt, ro, vro, a, mu):
    """
    Solves the universal Kepler equation using Newton's method.
    
    Parameters:
    dt : float
        Time since x = 0 (s).
    ro : float
        Radial position (km) when x = 0.
    vro : float
        Radial velocity (km/s) when x = 0.
    a : float
        Reciprocal of the semimajor axis (1/km).
        
    Returns:
    x : float
        The universal anomaly (km^0.5).
    """
    
    
    # Set an error tolerance and maximum number of iterations
    error = 1.e-8
    nMax = 1000
    
    # Starting value for x
    x = math.sqrt(mu) * abs(a) * dt
    
    # Iteration on Equation 3.65 until convergence or maximum iterations reached
    n = 0
    ratio = 1.0
    
    while abs(ratio) > error and n <= nMax:
        n += 1
        
        # Compute Stumpff functions and derivatives
        z = a * x**2
        C = stumpC(z)
        S = stumpS(z)
        
        # Calculate F(x) and its derivative dF/dx
        F = ro * vro / math.sqrt(mu) * x**2 * C + (1 - a * ro) * x**3 * S + ro * x - math.sqrt(mu) * dt
        dFdx = ro * vro / math.sqrt(mu) * x * (1 - a * x**2 * S) + (1 - a * ro) * x**2 * C + ro
        
        # Update x using Newton's method
        ratio = F / dFdx
        x = x - ratio
    
    # Report if maximum iterations were reached
    if n > nMax:
        print(f"\n **No. iterations of Kepler's equation = {n}")
        print(f" F/dFdx = {F/dFdx}")
    
    return x