import numpy as np

def q_from_dcm(Q):
    """
    Calculates the quaternion from the direction cosine matrix (DCM).

    Parameters:
    Q (numpy array): Direction cosine matrix (3x3).

    Returns:
    q (numpy array): Quaternion [q1, q2, q3, q4] where q4 is the scalar part.
    """
    K3 = np.array([
        [Q[0, 0] - Q[1, 1] - Q[2, 2], Q[1, 0] + Q[0, 1], Q[2, 0] + Q[0, 2], Q[1, 2] - Q[2, 1]],
        [Q[1, 0] + Q[0, 1], Q[1, 1] - Q[0, 0] - Q[2, 2], Q[2, 1] + Q[1, 2], Q[2, 0] - Q[0, 2]],
        [Q[2, 0] + Q[0, 2], Q[2, 1] + Q[1, 2], Q[2, 2] - Q[0, 0] - Q[1, 1], Q[0, 1] - Q[1, 0]],
        [Q[1, 2] - Q[2, 1], Q[2, 0] - Q[0, 2], Q[0, 1] - Q[1, 0], Q[0, 0] + Q[1, 1] + Q[2, 2]]
    ]) / 3.0

    eigval, eigvec = np.linalg.eig(K3)
    max_index = np.argmax(eigval)
    q = eigvec[:, max_index]

    if q[3] < 0:
        q = -q

    return q
