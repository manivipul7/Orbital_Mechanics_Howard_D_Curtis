import numpy as np
import matplotlib.pyplot as plt
from rkf45 import rkf45

# Constants and Input Data
days = 24 * 3600
G = 6.6742e-20  # km^3/kg/s^2
rmoon = 1737  # km
rearth = 6378  # km
r12 = 384400  # km
m1 = 5.974e24  # kg
m2 = 7.348e22  # kg
M = m1 + m2  # kg

pi_1 = m1 / M  # Moon's mass fraction
pi_2 = m2 / M  # Earth's mass fraction

mu1 = 398600  # Earth's gravitational parameter
mu2 = 4903.02  # Moon's gravitational parameter
mu = mu1 + mu2  # Earth-Moon system gravitational parameter

W = np.sqrt(mu / r12**3)  # Angular velocity of the moon around the Earth (rad/s)
x1 = -pi_2 * r12  # x-coordinate of the Earth relative to Earth-Moon barycenter
x2 = pi_1 * r12  # x-coordinate of the Moon relative to Earth-Moon barycenter

# Input data
d0 = 200  # Initial altitude of spacecraft (km)
phi = -90  # Initial polar azimuth coordinate (degrees)
v0 = 10.9148  # Initial speed of spacecraft relative to rotating Earth-Moon system (km/s)
gamma = 20  # Initial flight path angle (degrees)
t0 = 0  # Initial time (s)
tf = 3.16689 * days  # Final time (s)

# Initial conditions
r0 = rearth + d0  # Initial radial distance of spacecraft from Earth's center (km)
phi_rad = np.deg2rad(phi)
x = r0 * np.cos(phi_rad) + x1  # Initial x-coordinate in rotating Earth-Moon system (km)
y = r0 * np.sin(phi_rad)  # Initial y-coordinate in rotating Earth-Moon system (km)
vx = v0 * (np.sin(np.deg2rad(gamma)) * np.cos(phi_rad) - np.cos(np.deg2rad(gamma)) * np.sin(phi_rad))
vy = v0 * (np.sin(np.deg2rad(gamma)) * np.sin(phi_rad) + np.cos(np.deg2rad(gamma)) * np.cos(phi_rad))
f0 = np.array([x, y, vx, vy])  # Initial state vector

# Define the rates function
def rates(t, f):
    """
    Calculates the derivatives of the state vector f at time t.
    
    Parameters:
    t -- current time (s)
    f -- state vector [x, y, vx, vy] at time t
    
    Returns:
    dfdt -- derivative of state vector [vx, vy, ax, ay] at time t
    """
    x, y, vx, vy = f
    r1 = np.linalg.norm([x + pi_2 * r12, y])
    r2 = np.linalg.norm([x - pi_1 * r12, y])
    ax = 2 * W * vy + W**2 * x - mu1 * (x - x1) / r1**3 - mu2 * (x - x2) / r2**3
    ay = -2 * W * vx + W**2 * y - (mu1 / r1**3 + mu2 / r2**3) * y
    dfdt = np.array([vx, vy, ax, ay])
    return dfdt

# Output function
def output(t, f):
    """
    Echoes the input data, prints the results, and plots the trajectory.
    
    Parameters:
    t -- array of time points
    f -- array of state vectors [x, y, vx, vy] at each time point
    """
    xf, yf, vxf, vyf = f[-1,0], f[-1,1], f[-1,2], f[-1,3]
    df = np.sqrt((xf - x2)**2 + yf**2) - rmoon
    vf = np.sqrt(vxf**2 + vyf**2)
    
    print('\n\n––––––––––––––––––––––––––––––––––––––--––––––\n')
    print('\n Example 2.18: Lunar trajectory using the restricted')
    print('\n three-body equations.\n')
    print(f'\n Initial Earth altitude (km): {d0}')
    print(f'\n Initial angle between radial and earth-moon line (degrees): {phi}')
    print(f'\n Initial flight path angle (degrees): {gamma}')
    print(f'\n Flight time (days): {tf/days}')
    print(f'\n Final distance from the moon (km): {df}')
    print(f'\n Final relative speed (km/s): {vf}')
    print('\n\n––––––––––––––––––––––––––––––––––––––--––––––\n')
    
    # Print x and y values
    print("Printing x and y values:")
    for point in f:
        x, y, _, _ = point
        print(f"x: {x}, y: {y}")

    # Plot the trajectory and place filled circles representing the earth and moon
    plt.figure()
    plt.plot(f[:, 0], f[:, 1])
    plt.xlabel('x, km')
    plt.ylabel('y, km')
    plt.title('Lunar Trajectory')
    
    # Plot the earth (blue) and moon (green) to scale
    earth = circle(x1, 0, rearth)
    moon = circle(x2, 0, rmoon)
    plt.fill(earth[:, 0], earth[:, 1], 'b')
    plt.fill(moon[:, 0], moon[:, 1], 'g')
    
    # Set plot display parameters
    xmin, xmax = -20e3, 4e5
    ymin, ymax = -20e3, 1e5
    plt.axis([xmin, xmax, ymin, ymax])
    plt.axis('equal')
    plt.grid(True)
    plt.show()

# Helper function to create a circle for plotting
def circle(xc, yc, radius):
    """
    Creates coordinates of points spaced 0.1 degree apart around a circle.
    
    Parameters:
    xc, yc -- center coordinates of the circle
    radius -- radius of the circle
    
    Returns:
    xy -- array containing x and y coordinates of points on the circle
    """
    theta = np.linspace(0, 360, 3601)
    x = xc + radius * np.cos(np.deg2rad(theta))
    y = yc + radius * np.sin(np.deg2rad(theta))
    xy = np.column_stack([x, y])
    return xy

# Main execution
if __name__ == '__main__':
    # Perform the numerical integration using RKF45
    t, f = rkf45(rates, [t0, tf], f0)
    
    # Output the results
    output(t, f)
