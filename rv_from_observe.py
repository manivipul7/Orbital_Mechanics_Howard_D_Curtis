import numpy as np

def rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H, Re, f, wE):
    # Constants
    deg = np.pi / 180
    omega = np.array([0, 0, wE])

    # Convert angular quantities from degrees to radians
    A = A * deg
    Adot = Adot * deg
    a = a * deg
    adot = adot * deg
    theta = theta * deg
    phi = phi * deg

    # Equation 5.56
    R = np.array([
        (Re / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.cos(phi) * np.cos(theta),
        (Re / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.cos(phi) * np.sin(theta),
        (Re * (1 - f)**2 / np.sqrt(1 - (2 * f - f**2) * np.sin(phi)**2) + H) * np.sin(phi)
    ])

    # Equation 5.66
    Rdot = np.cross(omega, R)

    # Equation 5.83a
    dec = np.arcsin(np.cos(phi) * np.cos(A) * np.cos(a) + np.sin(phi) * np.sin(a))

    # Equation 5.83b
    h = np.arccos((np.cos(phi) * np.sin(a) - np.sin(phi) * np.cos(A) * np.cos(a)) / np.cos(dec))
    if 0 < A < np.pi:
        h = 2 * np.pi - h

    # Equation 5.83c
    RA = theta - h

    # Equation 5.57
    Rho = np.array([np.cos(RA) * np.cos(dec), np.sin(RA) * np.cos(dec), np.sin(dec)])

    # Equation 5.63
    r = R + rho * Rho

    # Equation 5.84
    decdot = (-Adot * np.cos(phi) * np.sin(A) * np.cos(a) + adot * (np.sin(phi) * np.cos(a) - np.cos(phi) * np.cos(A) * np.sin(a))) / np.cos(dec)

    # Equation 5.85
    RAdot = wE + (Adot * np.cos(A) * np.cos(a) - adot * np.sin(A) * np.sin(a) + decdot * np.sin(A) * np.cos(a) * np.tan(dec)) / (np.cos(phi) * np.sin(a) - np.sin(phi) * np.cos(A) * np.cos(a))

    # Equations 5.69 and 5.72
    Rhodot = np.array([
        -RAdot * np.sin(RA) * np.cos(dec) - decdot * np.cos(RA) * np.sin(dec),
        RAdot * np.cos(RA) * np.cos(dec) - decdot * np.sin(RA) * np.sin(dec),
        decdot * np.cos(dec)
    ])

    # Equation 5.64
    v = Rdot + rhodot * Rho + rho * Rhodot

    return r, v


