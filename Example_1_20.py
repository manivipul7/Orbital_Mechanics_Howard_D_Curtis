import numpy as np
import matplotlib.pyplot as plt
from rkf45 import rkf45

def rates(t, f):
    """
    This function calculates first and second time derivatives of x
    governed by the equation of two-body rectilinear motion.
    x'' + mu/x^2 = 0
    """
    x = f[0]
    Dx = f[1]
    D2x = -mu / x**2
    return np.array([Dx, D2x])

def plot_results(t, f):
    minutes = 60  # Conversion from seconds to minutes
    
    # Position vs time
    plt.subplot(2, 1, 1)
    plt.plot(t / minutes, f[:, 0], '-ok')
    plt.xlabel('time, minutes')
    plt.ylabel('position, km')
    plt.grid(True)
    plt.axis([0, np.max(t)/minutes, 5000, 15000])
    
    # Velocity vs time
    plt.subplot(2, 1, 2)
    plt.plot(t / minutes, f[:, 1], '-ok')
    plt.xlabel('time, minutes')
    plt.ylabel('velocity, km/s')
    plt.grid(True)
    plt.axis([0, np.max(t)/minutes, -10, 10])
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # Define constants
    mu = 398600  # km^3/s^2
    minutes = 60  # Conversion from minutes to seconds
    
    # Initial conditions
    x0 = 6500  # Initial position (km)
    v0 = 7.8   # Initial velocity (km/s)
    y0 = [x0, v0]  # Initial state vector
    t0 = 0     # Initial time (s)
    tf = 70 * minutes  # Final time (s)
    tspan = [t0, tf]
    
    # Solve the differential equation using RKF45
    t, f = rkf45(rates, tspan, y0)
    
    # Plot the results
    plot_results(t, f)
