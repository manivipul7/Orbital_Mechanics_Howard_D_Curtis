import numpy as np
from atan2d_0_360 import atan2d_0_360

def dcm_to_ypr(Q):
    """
    This function finds the angles of the yaw-pitch-roll sequence
    R1(gamma)*R2(beta)*R3(alpha) from the direction cosine matrix.
    
    Parameters:
    Q -- direction cosine matrix (3x3 numpy array)
    
    Returns:
    yaw -- yaw angle (degrees)
    pitch -- pitch angle (degrees)
    roll -- roll angle (degrees)
    """
    yaw = atan2d_0_360(Q[0, 1], Q[0, 0])
    pitch = np.degrees(np.arcsin(-Q[0, 2]))
    roll = atan2d_0_360(Q[1, 2], Q[2, 2])
    
    return yaw, pitch, roll