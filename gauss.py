from kepler_U import kepler_U
from f_and_g import f_and_g
import numpy as np
from scipy.optimize import fsolve

def gauss(Rho1, Rho2, Rho3, R1, R2, R3, t1, t2, t3, mu):
    # Compute time intervals between observations
    tau1 = t1 - t2
    tau3 = t3 - t2
    tau = tau3 - tau1

    # Compute cross products among the direction cosine vectors
    p1 = np.cross(Rho2, Rho3)
    p2 = np.cross(Rho1, Rho3)
    p3 = np.cross(Rho1, Rho2)

    # Compute Do (dot product of Rho1 and p1)
    Do = np.dot(Rho1, p1)

    # Compute the matrix D of the nine scalar triple products
    D = np.array([
        [np.dot(R1, p1), np.dot(R1, p2), np.dot(R1, p3)],
        [np.dot(R2, p1), np.dot(R2, p2), np.dot(R2, p3)],
        [np.dot(R3, p1), np.dot(R3, p2), np.dot(R3, p3)]
    ])

    # Compute E (dot product of R2 and Rho2)
    E = np.dot(R2, Rho2)

    # Compute A and B
    A = 1 / Do * (-D[0, 1] * tau3 / tau + D[1, 1] + D[2, 1] * tau1 / tau)
    B = 1 / (6 * Do) * (D[0, 1] * (tau3**2 - tau**2) * tau3 / tau + D[2, 1] * (tau**2 - tau1**2) * tau1 / tau)

    # Compute coefficients of the 8th order polynomial
    a_coeff = -(A**2 + 2 * A * E + np.linalg.norm(R2)**2)
    b_coeff = -2 * mu * B * (A + E)
    c_coeff = -(mu * B)**2

    # Solve for the positive real root
    Roots = np.roots([1, 0, a_coeff, 0, 0, b_coeff, 0, 0, c_coeff])
    x = posroot(Roots)

    # Compute f1, f3, g1, g3
    f1 = 1 - 1/2 * mu * tau1**2 / x**3
    f3 = 1 - 1/2 * mu * tau3**2 / x**3
    g1 = tau1 - 1/6 * mu * (tau1 / x)**3
    g3 = tau3 - 1/6 * mu * (tau3 / x)**3

    # Compute rho2, rho1, rho3
    rho2 = A + mu * B / x**3
    rho1 = 1 / Do * ((6 * (D[2, 0] * tau1 / tau3 + D[1, 0] * tau / tau3) * x**3 + mu * D[2, 0] * (tau**2 - tau1**2) * tau1 / tau3) / (6 * x**3 + mu * (tau**2 - tau3**2)) - D[0, 0])
    rho3 = 1 / Do * ((6 * (D[0, 2] * tau3 / tau1 - D[1, 2] * tau / tau1) * x**3 + mu * D[0, 2] * (tau**2 - tau3**2) * tau3 / tau1) / (6 * x**3 + mu * (tau**2 - tau1**2)) - D[2, 2])

    # Compute position vectors
    r1 = R1 + rho1 * Rho1
    r2 = R2 + rho2 * Rho2
    r3 = R3 + rho3 * Rho3

    # Compute velocity vector
    v2 = (-f3 * r1 + f1 * r3) / (f1 * g3 - f3 * g1)

    # Save the initial estimates
    r_old = r2.copy()
    v_old = v2.copy()

    # Iterative improvement loop
    rho1_old, rho2_old, rho3_old = rho1, rho2, rho3
    tol = 1e-8
    nmax = 1000
    n = 0
    diff1, diff2, diff3 = 1, 1, 1

    while ((diff1 > tol) or (diff2 > tol) or (diff3 > tol)) and (n < nmax):
        n += 1

        # Compute quantities required by universal Kepler's equation
        ro = np.linalg.norm(r2)
        vo = np.linalg.norm(v2)
        vro = np.dot(v2, r2) / ro
        a = 2 / ro - vo**2 / mu

        # Solve universal Kepler's equation
        x1 = kepler_U(tau1, ro, vro, a, mu)
        x3 = kepler_U(tau3, ro, vro, a, mu)

        # Calculate the Lagrange f and g coefficients
        f1_new, g1_new = f_and_g(x1, tau1, ro, a, mu)
        f3_new, g3_new = f_and_g(x3, tau3, ro, a, mu)

        # Update f and g functions
        f1 = (f1 + f1_new) / 2
        f3 = (f3 + f3_new) / 2
        g1 = (g1 + g1_new) / 2
        g3 = (g3 + g3_new) / 2

        # Update rho1, rho2, rho3
        c1 = g3 / (f1 * g3 - f3 * g1)
        c3 = -g1 / (f1 * g3 - f3 * g1)

        rho1 = 1 / Do * (-D[0, 0] + 1 / c1 * D[1, 0] - c3 / c1 * D[2, 0])
        rho2 = 1 / Do * (-c1 * D[0, 1] + D[1, 1] - c3 * D[2, 1])
        rho3 = 1 / Do * (-c1 / c3 * D[0, 2] + 1 / c3 * D[1, 2] - D[2, 2])

        # Update position vectors
        r1 = R1 + rho1 * Rho1
        r2 = R2 + rho2 * Rho2
        r3 = R3 + rho3 * Rho3

        # Update velocity vector
        v2 = (-f3 * r1 + f1 * r3) / (f1 * g3 - f3 * g1)

        # Calculate differences
        diff1 = abs(rho1 - rho1_old)
        diff2 = abs(rho2 - rho2_old)
        diff3 = abs(rho3 - rho3_old)

        # Update slant ranges
        rho1_old, rho2_old, rho3_old = rho1, rho2, rho3

    # Print the number of iterations
    print(f"\n** Number of Gauss improvement iterations = {n}\n")
    if n >= nmax:
        print(f"\n** Number of iterations exceeds {nmax}\n")

    # Return the state vector for the central observation
    return r2, v2, r_old, v_old

def posroot(Roots):
    posroots = Roots[np.where((Roots > 0) & (np.imag(Roots) == 0))]
    if len(posroots) == 0:
        raise ValueError("** There are no positive roots.")
    elif len(posroots) == 1:
        return posroots[0].real
    else:
        for i, root in enumerate(posroots):
            print(f"root #{i + 1} = {root.real}")
        nchoice = int(input("Use root #? ")) - 1
        if 0 <= nchoice < len(posroots):
            return posroots[nchoice].real
        else:
            raise ValueError("Invalid choice.")

# Example implementations of kepler_U and f_and_g functions are required
# Define your kepler_U and f_and_g functions here.
# Here are placeholders for these functions:
