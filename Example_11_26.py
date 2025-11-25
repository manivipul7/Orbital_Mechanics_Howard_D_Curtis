import numpy as np
import matplotlib.pyplot as plt
from q_from_dcm import q_from_dcm
from rkf45 import rkf45
from dcm_from_q import dcm_from_q
from dcm_to_euler import dcm_to_euler

def Example_11_26():
    """
    Numerically integrates Euler’s equations of motion for a spinning top.
    Uses quaternions to obtain the time history of the top’s orientation.
    """
    global m, g, d, A, B, C, ws0, wp0, wn0, theta
    # Constants and initial data
    g = 9.807  # Acceleration of gravity (m/s^2)
    m = 0.5    # Mass (kg)
    d = 0.05   # Distance of center of mass from pivot point (m)
    A = 12e-4  # Moment of inertia about body x (kg-m^2)
    B = 12e-4  # Moment of inertia about body y (kg-m^2)
    C = 4.5e-4 # Moment of inertia about body z (kg-m^2)
    ws0 = 1000 * 2 * np.pi / 60  # Spin rate (rad/s)
    wp0 = 0    # Initial precession rate (rad/s)
    wn0 = 0    # Initial nutation rate (rad/s)
    theta = 60 # Initial nutation angle (deg)

    z = np.array([0, -np.sin(np.deg2rad(theta)), np.cos(np.deg2rad(theta))])
    p = np.array([1, 0, 0])

    # Initial direction cosine matrix
    y = np.cross(z, p)
    x = np.cross(y, z)
    i = x / np.linalg.norm(x)
    j = y / np.linalg.norm(y)
    k = z / np.linalg.norm(z)
    QXx = np.array([i, j, k])

    # Initial quaternion
    q0 = q_from_dcm(QXx)

    phi0, theta0, psi0 = dcm_to_euler(QXx)

    # Initial body-frame angular velocity
    w0 = np.array([wp0 * np.sin(np.deg2rad(theta0)) * np.sin(np.deg2rad(psi0)) + wn0 * np.cos(np.deg2rad(psi0)),
                   wp0 * np.sin(np.deg2rad(theta0)) * np.cos(np.deg2rad(psi0)) - wn0 * np.sin(np.deg2rad(psi0)),
                   ws0 + wp0 * np.cos(np.deg2rad(theta0))])

    # Initial conditions vector
    f0 = np.concatenate((q0, w0))

    # Integration time
    t0 = 0
    tf = 1.153  # Final time for 360 degrees of precession

    # RKF4(5) numerical ODE solver
    t, f = rkf45(rates, [t0, tf], f0)

    # Solutions for quaternion and angular velocities
    q = f[:, 0:4]
    wx = f[:, 4]
    wy = f[:, 5]
    wz = f[:, 6]

    # Arrays to store results
    prec = np.zeros_like(t)
    nut = np.zeros_like(t)
    spin = np.zeros_like(t)
    wp = np.zeros_like(t)
    wn = np.zeros_like(t)
    ws = np.zeros_like(t)

    # Obtain direction cosine matrix, Euler angles, and Euler angle rates
    for m in range(len(t)):
        QXx = dcm_from_q(q[m])
        prec[m], nut[m], spin[m] = dcm_to_euler(QXx)
        wp[m] = (wx[m] * np.sin(np.deg2rad(spin[m])) + wy[m] * np.cos(np.deg2rad(spin[m]))) / np.sin(np.deg2rad(nut[m]))
        wn[m] = wx[m] * np.cos(np.deg2rad(spin[m])) - wy[m] * np.sin(np.deg2rad(spin[m]))
        ws[m] = -wp[m] * np.cos(np.deg2rad(nut[m])) + wz[m]

    # Plotting results
    plotit(t, prec, wp, nut, wn, spin, ws)

def rates(t, f):
    """
    Time derivative function for the spinning top Euler's equations of motion.

    Parameters:
    t (float): Time parameter (s).
    f (numpy array): State vector containing quaternion and angular velocities.

    Returns:
    dfdt (numpy array): Time derivatives of the state vector.
    """
    q = f[0:4]    # Quaternion components
    wx = f[4]     # Angular velocity x-component
    wy = f[5]     # Angular velocity y-component
    wz = f[6]     # Angular velocity z-component

    q /= np.linalg.norm(q)  # Normalize quaternion

    # Angular velocity vector
    w = np.array([wx, wy, wz])

    # Direction cosine matrix from quaternion
    Q = dcm_from_q(q)

    # Moment components about the pivot point
    M = np.matmul(Q, np.array([-m * g * d * Q[2, 1], 
                                m * g * d * Q[2, 0], 
                                                0]))

    # Skew-symmetric matrix of angular velocities
    Omega = np.array([[0, wz, -wy, wx],
                      [-wz, 0, wx, wy],
                      [wy, -wx, 0, wz],
                      [-wx, -wy, -wz, 0]])

    # Time derivative of quaternion
    q_dot = 0.5 * np.matmul(Omega, q)

    # Euler's equations
    wx_dot = M[0] / A - (C - B) * wy * wz / A
    wy_dot = M[1] / B - (A - C) * wz * wx / B
    wz_dot = M[2] / C - (B - A) * wx * wy / C

    # Time derivatives vector
    dfdt = np.concatenate((q_dot, [wx_dot, wy_dot, wz_dot]))

    return dfdt

def plotit(t, prec, wp, nut, wn, spin, ws):
    """
    Function to plot Euler angles and their rates.

    Parameters:
    t (numpy array): Time array (s).
    prec (numpy array): Precession angle (deg).
    wp (numpy array): Precession rate (rpm).
    nut (numpy array): Nutation angle (deg).
    wn (numpy array): Nutation rate (deg/s).
    spin (numpy array): Spin angle (deg).
    ws (numpy array): Spin rate (rpm).
    """
    plt.figure(figsize=(10, 12))
    plt.subplot(321)
    plt.plot(t, prec)
    plt.xlabel('Time (s)')
    plt.ylabel('Precession angle (deg)')
    plt.grid(True)

    plt.subplot(322)
    plt.plot(t, wp * 60 / 2 / np.pi)
    plt.xlabel('Time (s)')
    plt.ylabel('Precession rate (rpm)')
    plt.grid(True)

    plt.subplot(323)
    plt.plot(t, nut)
    plt.xlabel('Time (s)')
    plt.ylabel('Nutation angle (deg)')
    plt.grid(True)

    plt.subplot(324)
    plt.plot(t, wn * 180 / np.pi)
    plt.xlabel('Time (s)')
    plt.ylabel('Nutation rate (deg/s)')
    plt.grid(True)

    plt.subplot(325)
    plt.plot(t, spin)
    plt.xlabel('Time (s)')
    plt.ylabel('Spin angle (deg)')
    plt.grid(True)

    plt.subplot(326)
    plt.plot(t, ws * 60 / 2 / np.pi)
    plt.xlabel('Time (s)')
    plt.ylabel('Spin rate (rpm)')
    plt.grid(True)

    plt.tight_layout()
    plt.show()

Example_11_26()
