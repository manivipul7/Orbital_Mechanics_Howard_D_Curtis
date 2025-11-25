from stumpC import stumpC
import numpy as np
from stumpS import stumpS

def lambert(R1, R2, t, string='pro'):
    global mu
    mu = 398600  # Gravitational parameter (km^3/s^2)

    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    c12 = np.cross(R1, R2)
    theta = np.arccos(np.dot(R1, R2) / (r1 * r2))

    if string == 'pro':
        if c12[2] <= 0:
            theta = 2 * np.pi - theta
    elif string == 'retro':
        if c12[2] >= 0:
            theta = 2 * np.pi - theta

    A = np.sin(theta) * np.sqrt(r1 * r2 / (1 - np.cos(theta)))

    z = 0.1
    while F(z, t, A, r1, r2, theta) < 0:
        z += 0.1

    tol = 1e-8
    nmax = 5000
    ratio = 1
    n = 0
    while abs(ratio) > tol and n <= nmax:
        n += 1
        Fz = F(z, t, A, r1, r2, theta)
        dFz_dz = dFdz(z, r1, r2, theta, A)
        ratio = Fz / dFz_dz
        z -= ratio

    if n >= nmax:
        print('Number of iterations exceeds maximum.')

    f = 1 - y(z, r1, r2, theta, A) / r1
    g = A * np.sqrt(y(z, r1, r2, theta, A) / mu)
    gdot = 1 - y(z, r1, r2, theta, A) / r2

    V1 = (R2 - f * R1) / g
    V2 = (gdot * R2 - R1) / g

    return V1, V2

def y(z, r1, r2, theta, A):
    global mu
    return r1 + r2 + A * ((z * stumpS(z) - 1) / np.sqrt(stumpC(z)))

def F(z, t, A, r1, r2, theta):
    global mu
    return (y(z, r1, r2, theta, A) / stumpC(z))**(3 / 2) * stumpS(z) + A * np.sqrt(y(z, r1, r2, theta, A)) - np.sqrt(mu) * t

def dFdz(z, r1, r2, theta, A):
    global mu
    if z == 0:
        return np.sqrt(2) / 40 * y(0, r1, r2, theta, A)**(3 / 2) + A / 8 * (np.sqrt(y(0, r1, r2, theta, A)) + A * np.sqrt(1 / (2 * y(0, r1, r2, theta, A))))
    else:
        return (y(z, r1, r2, theta, A) / stumpC(z))**(3 / 2) * (1 / (2 * z) * (stumpC(z) - 3 * stumpS(z) / (2 * stumpC(z))) + 
                (3 * stumpS(z)**2) / 4 / stumpC(z)) + A / 8 * (3 * stumpS(z) / stumpC(z) * np.sqrt(y(z, r1, r2, theta, A)) + 
                A * np.sqrt(stumpC(z) / y(z, r1, r2, theta, A)))
