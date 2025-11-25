import numpy as np

def quat_rotate(q, v):
    """
    Rotates a vector `v` by a unit quaternion `q`.

    Parameters:
    q (numpy array): Unit quaternion [q1, q2, q3, q4] where q4 is the scalar part.
    v (numpy array): Vector to be rotated [vx, vy, vz].

    Returns:
    r (numpy array): Rotated vector [rx, ry, rz].
    """
    qinv = np.array([q[0], -q[1], -q[2], -q[3]])  # Quaternion inverse

    # Convert v to a pure quaternion [0 v]
    V = np.array([0, v[0], v[1], v[2]])

    # Quaternion multiplication q * V * qinv
    R = quatmultiply(quatmultiply(q, V), qinv)

    # Extract the rotated vector [R(2), R(3), R(4)]
    r = R[1:]

    return r

def quatmultiply(q1, q2):
    """
    Multiplies two quaternions.

    Parameters:
    q1 (numpy array): Quaternion [q1, q2, q3, q4] where q4 is the scalar part.
    q2 (numpy array): Quaternion [q1, q2, q3, q4] where q4 is the scalar part.

    Returns:
    q (numpy array): Product of quaternions [q1*q2, q1*q3, q1*q4, q4*q2].
    """
    q = np.zeros(4)
    q[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3]
    q[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2]
    q[2] = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1]
    q[3] = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0]
    return q
