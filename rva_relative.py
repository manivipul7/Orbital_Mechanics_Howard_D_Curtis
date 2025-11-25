import numpy as np

def rva_relative(rA, vA, rB, vB, mu):
    # Calculate the vector hA
    hA = np.cross(rA, vA)

    # Calculate the unit vectors i, j, and k
    i = rA / np.linalg.norm(rA)
    k = hA / np.linalg.norm(hA)
    j = np.cross(k, i)

    # Calculate the transformation matrix QXx
    QXx = np.array([i, j, k])

    # Calculate Omega
    Omega = hA / np.linalg.norm(rA)**2

    # Calculate Omega_dot
    Omega_dot = -2 * np.dot(rA, vA) / np.linalg.norm(rA)**2 * Omega

    # Calculate the accelerations aA and aB
    aA = -mu * rA / np.linalg.norm(rA)**3
    aB = -mu * rB / np.linalg.norm(rB)**3

    # Calculate r_rel, v_rel, a_rel
    r_rel = rB - rA
    v_rel = vB - vA - np.cross(Omega, r_rel)
    a_rel = aB - aA - np.cross(Omega_dot, r_rel) \
            - np.cross(Omega, np.cross(Omega, r_rel)) \
            - 2 * np.cross(Omega, v_rel)

    # Calculate r_rel_x, v_rel_x, a_rel_x
    r_rel_x = np.dot(QXx, r_rel)
    v_rel_x = np.dot(QXx, v_rel)
    a_rel_x = np.dot(QXx, a_rel)

    return r_rel_x, v_rel_x, a_rel_x
