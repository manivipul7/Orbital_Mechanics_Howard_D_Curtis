import numpy as np
from my_bisect import bisect

def fun(z, p):
    """
    This function evaluates the function in Equation 2.204.

    Parameters:
    z -- the dimensionless x-coordinate
    p -- ratio of moon mass to total mass

    Returns:
    f -- the value of the function
    """
    f = (1 - p) * (z + p) / np.abs(z + p)**3 + p * (z + p - 1) / np.abs(z + p - 1)**3 - z
    return f

def output(x, p, r12):
    """
    This function prints out the x coordinates of L1, L2, and L3 relative to the center of mass.

    Parameters:
    x -- vector containing the three computed roots
    p -- ratio of moon mass to total mass
    r12 -- distance from the earth to the moon (km)
    """
    print('\n\n––––––––––––––––––––––––––––––––––––––--––––––\n')
    print('\n For\n')
    print(f'\n m1 = {m1:.3e} kg')
    print(f'\n m2 = {m2:.3e} kg')
    print(f'\n r12 = {r12:.3e} km\n')
    print('\n the 3 colinear Lagrange points (the roots of\n')
    print(' Equation 2.204) are:\n')
    print(f'\n L3: x = {x[0]*r12:.10g} km (f(x3) = {fun(x[0], p):.10g})')
    print(f'\n L1: x = {x[1]*r12:.10g} km (f(x1) = {fun(x[1], p):.10g})')
    print(f'\n L2: x = {x[2]*r12:.10g} km (f(x2) = {fun(x[2], p):.10g})')
    print('\n\n––––––––––––––––––––––––––––––––––––––--––––––\n')

if __name__ == "__main__":
    # Input data
    m1 = 5.974e24  # mass of the earth (kg)
    m2 = 7.348e22  # mass of the moon (kg)
    r12 = 3.844e5  # distance from the earth to the moon (km)
    xl = np.array([-1.1, 0.5, 1.0])  # low-side estimates of the three roots
    xu = np.array([-0.9, 1.0, 1.5])  # high-side estimates of the three roots

    # Compute ratio of moon mass to total mass
    p = m2 / (m1 + m2)

    # Compute the roots using the bisection method
    x = np.zeros(3)
    for i in range(3):
        x[i] = bisect(lambda z: fun(z, p), xl[i], xu[i])

    # Output the results
    output(x, p, r12)
