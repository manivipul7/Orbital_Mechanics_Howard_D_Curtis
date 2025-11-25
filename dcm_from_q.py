import numpy as np

def dcm_from_q(q):
    """
    Calculates the direction cosine matrix (DCM) from a quaternion.

    Parameters:
    q (numpy array): Quaternion [q1, q2, q3, q4] where q4 is the scalar part.

    Returns:
    Q (numpy array): Direction cosine matrix (3x3).
    """
    q1, q2, q3, q4 = q[0], q[1], q[2], q[3]

    Q = np.array([
        [q1**2 - q2**2 - q3**2 + q4**2, 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4)],
        [2*(q1*q2 - q3*q4), -q1**2 + q2**2 - q3**2 + q4**2, 2*(q2*q3 + q1*q4)],
        [2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), -q1**2 - q2**2 + q3**2 + q4**2]
    ])

    return Q
