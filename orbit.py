import numpy as np
import matplotlib.pyplot as plt
from rkf45 import rkf45
from mpl_toolkits.mplot3d import Axes3D

def rates(t, f, mu):
    """
    This function calculates the acceleration vector using the two-body problem equations.
    """
    x, y, z = f[0], f[1], f[2]
    vx, vy, vz = f[3], f[4], f[5]
    r = np.linalg.norm([x, y, z])
    ax = -mu * x / r**3
    ay = -mu * y / r**3
    az = -mu * z / r**3
    dydt = np.array([vx, vy, vz, ax, ay, az])
    return dydt

def light_gray():
    """
    This function creates a color map for displaying the planet as light gray with a black equator.
    """
    r = 0.8
    g = r
    b = r
    return np.array([[r, g, b], [0, 0, 0], [r, g, b]])

def output(t, y, r0, v0, R):
    """
    This function computes the maximum and minimum radii, the times they occur, and the speed at those times.
    It prints those results and plots the orbit.
    """
    r = np.linalg.norm(y[:, :3], axis=1)
    imax = np.argmax(r)
    imin = np.argmin(r)
    rmax = r[imax]
    rmin = r[imin]
    v_at_rmax = np.linalg.norm(y[imax, 3:])
    v_at_rmin = np.linalg.norm(y[imin, 3:])

    # Output to the command window:
    print("\n\n––––––––––––––––––––––––––––––––––––––--––––––––––––––––––––––\n")
    print("\n Earth Orbit\n")
    print(f" The initial position is [{r0[0]}, {r0[1]}, {r0[2]}] (km).")
    print(f" Magnitude = {np.linalg.norm(r0)} km\n")
    print(f" The initial velocity is [{v0[0]}, {v0[1]}, {v0[2]}] (km/s).")
    print(f" Magnitude = {np.linalg.norm(v0)} km/s\n")
    print(f" Initial time = {0} h.\n Final time = {tf / 3600} h.\n")
    print(f" The minimum altitude is {rmin - R} km at time = {t[imin] / 3600} h.")
    print(f" The speed at that point is {v_at_rmin} km/s.\n")
    print(f" The maximum altitude is {rmax - R} km at time = {t[imax] / 3600} h.")
    print(f" The speed at that point is {v_at_rmax} km/s\n")
    print("––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––\n\n")

    # Plot the results:
    # Draw the planet
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    xx = R * np.cos(u) * np.sin(v)
    yy = R * np.sin(u) * np.sin(v)
    zz = R * np.cos(v)
    ax.plot_surface(xx, yy, zz, color='lightgray', rstride=5, cstride=5, alpha=0.5, linewidth=0)

    # Draw and label the X, Y, and Z axes
    ax.plot([0, 2*R], [0, 0], [0, 0], color='k')
    ax.text(2*R, 0, 0, 'X')
    ax.plot([0, 0], [0, 2*R], [0, 0], color='k')
    ax.text(0, 2*R, 0, 'Y')
    ax.plot([0, 0], [0, 0], [0, 2*R], color='k')
    ax.text(0, 0, 2*R, 'Z')

    # Plot the orbit, draw a radial to the starting point, and label the starting point (o) and the final point (f)
    ax.plot(y[:, 0], y[:, 1], y[:, 2], 'k')
    ax.plot([0, r0[0]], [0, r0[1]], [0, r0[2]], 'b--')
    ax.text(y[0, 0], y[0, 1], y[0, 2], 'o')
    ax.text(y[-1, 0], y[-1, 1], y[-1, 2], 'f')

    # Select a view direction (a vector directed outward from the origin)
    ax.view_init(elev=24, azim=45)

    # Specify some properties of the graph
    ax.grid(True)
    ax.set_aspect('equal')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    plt.show()

if __name__ == '__main__':
    hours = 3600
    G = 6.6742e-20

    # Input data:
    m1 = 5.974e24  # Earth mass (kg)
    R = 6378  # Earth radius (km)
    m2 = 1000  # Spacecraft mass (kg)
    r0 = np.array([8000, 0, 6000])  # Initial position vector (km)
    v0 = np.array([0, 7, 0])  # Initial velocity vector (km/s)
    t0 = 0  # Initial time (s)
    tf = 4 * hours  # Final time (s)

    # Gravitational parameter
    mu = G * (m1 + m2)

    # Initial state vector
    y0 = np.concatenate((r0, v0))

    # Solve the differential equations using RKF45
    tspan = [t0, tf]
    t, y = rkf45(lambda t, f: rates(t, f, mu), tspan, y0)

    # Plot the results
    output(t, y, r0, v0, R)
