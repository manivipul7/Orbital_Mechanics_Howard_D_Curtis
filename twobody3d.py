import numpy as np
import matplotlib.pyplot as plt
from rkf45 import rkf45

def rates(t, y, G, m1, m2):
    """
    This function calculates the accelerations in the two-body problem.
    """
    R1 = y[0:3]
    R2 = y[3:6]
    V1 = y[6:9]
    V2 = y[9:12]
    r = np.linalg.norm(R2 - R1)
    A1 = G * m2 * (R2 - R1) / r**3
    A2 = G * m1 * (R1 - R2) / r**3
    dydt = np.concatenate((V1, V2, A1, A2))
    return dydt

def common_axis_settings(ax):
    """
    This function establishes axis properties common to the several plots.
    """
    ax.text(0, 0, 0, 'o')
    ax.grid(True)
    ax.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.view_init(elev=24, azim=45)

def plot_results(t, y, m1, m2):
    """
    This function calculates the trajectory of the center of mass and plots the results.
    """
    X1, Y1, Z1 = y[:, 0], y[:, 1], y[:, 2]
    X2, Y2, Z2 = y[:, 3], y[:, 4], y[:, 5]

    # Locate the center of mass at each time step:
    XG = (m1 * X1 + m2 * X2) / (m1 + m2)
    YG = (m1 * Y1 + m2 * Y2) / (m1 + m2)
    ZG = (m1 * Z1 + m2 * Z2) / (m1 + m2)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.set_title('Figure 2.3: Motion relative to the inertial frame')
    ax1.plot(X1, Y1, Z1, '-r', label='m1')
    ax1.plot(X2, Y2, Z2, '-g', label='m2')
    ax1.plot(XG, YG, ZG, '-b', label='Center of Mass')
    common_axis_settings(ax1)
    ax1.legend()

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.set_title('Figure 2.4a: Motion of m2 and G relative to m1')
    ax2.plot(X2 - X1, Y2 - Y1, Z2 - Z1, '-g', label='m2 relative to m1')
    ax2.plot(XG - X1, YG - Y1, ZG - Z1, '-b', label='Center of Mass relative to m1')
    common_axis_settings(ax2)
    ax2.legend()

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111, projection='3d')
    ax3.set_title('Figure 2.4b: Motion of m1 and m2 relative to G')
    ax3.plot(X1 - XG, Y1 - YG, Z1 - ZG, '-r', label='m1 relative to Center of Mass')
    ax3.plot(X2 - XG, Y2 - YG, Z2 - ZG, '-g', label='m2 relative to Center of Mass')
    common_axis_settings(ax3)
    ax3.legend()

    plt.show()

if __name__ == '__main__':
    # Define constants
    G = 6.67259e-20  # km^3/kg/s^2
    m1 = 1.e26  # kg
    m2 = 1.e26  # kg
    t0 = 0  # Initial time (s)
    tf = 480  # Final time (s)

    # Initial conditions
    R1_0 = np.array([0, 0, 0])  # Initial position of m1 (km)
    R2_0 = np.array([3000, 0, 0])  # Initial position of m2 (km)
    V1_0 = np.array([10, 20, 30])  # Initial velocity of m1 (km/s)
    V2_0 = np.array([0, 40, 0])  # Initial velocity of m2 (km/s)
    y0 = np.concatenate((R1_0, R2_0, V1_0, V2_0))

    # Define the time span for the integration
    tspan = [t0, tf]

    # Solve the differential equations using RKF45
    t, y = rkf45(lambda t, y: rates(t, y, G, m1, m2), tspan, y0)

    # Plot the results
    plot_results(t, y, m1, m2)
