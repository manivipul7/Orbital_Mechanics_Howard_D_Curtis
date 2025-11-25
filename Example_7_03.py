import numpy as np
import matplotlib.pyplot as plt
from sv_from_coe import sv_from_coe
from rv_from_r0v0 import rv_from_r0v0
from rkf45 import rkf45

# Define global constants
mu = 398600  # gravitational parameter (km^3/s^2)
RE = 6378    # Earth's radius (km)

def Example_7_03():
    global mu, RE
    
    # Input data: Prescribed initial orbital parameters of target A
    rp = RE + 300
    e = 0.1
    i = 0
    RA = 0
    omega = 0
    theta = 0
    
    # Additional computed parameters
    ra = rp * (1 + e) / (1 - e)
    h = np.sqrt(2 * mu * rp * ra / (ra + rp))
    a = (rp + ra) / 2
    T = 2 * np.pi / np.sqrt(mu) * a**1.5
    n = 2 * np.pi / T
    
    # Prescribed initial state vector of chaser B in the co-moving frame
    dr0 = np.array([-1, 0, 0])
    dv0 = np.array([0, -2 * n * dr0[0], 0])
    t0 = 0
    tf = 5 * T
    
    # Calculate the target’s initial state vector using Algorithm 4.5
    R0, V0 = sv_from_coe([h, e, RA, i, omega, theta], mu)
    
    # Initial state vector of B’s orbit relative to A
    y0 = np.concatenate((dr0, dv0))
    
    # Integrate Equations 7.34 using RK4 method (similar to rkf45 in MATLAB)
    t, y = rkf45(lambda t, y: rates(t, y, R0, V0), [t0, tf], y0)
    
    # Plot the results
    plotit(y)

def rates(t, f, R0, V0):
    # Function to compute the components of f(t, y) in Equation 7.36
    # This function should return dydt
    # Update R, V using Algorithm 3.4
    R, V = rv_from_r0v0(R0, V0, t)
    
    X, Y, Z = R
    VX, VY, VZ = V
    R_ = np.linalg.norm([X, Y, Z])
    RdotV = np.dot([X, Y, Z], [VX, VY, VZ])
    h = np.linalg.norm(np.cross([X, Y, Z], [VX, VY, VZ]))
    
    dx, dy, dz = f[:3]
    dvx, dvy, dvz = f[3:]
    
    dax = (2 * mu / R_**3 + h**2 / R_**4) * dx - 2 * RdotV / R_**4 * h * dy + 2 * h / R_**2 * dvy
    day = -(mu / R_**3 - h**2 / R_**4) * dy + 2 * RdotV / R_**4 * h * dx - 2 * h / R_**2 * dvx
    daz = -mu / R_**3 * dz
    
    return np.array([dvx, dvy, dvz, dax, day, daz])

def plotit(y):
    # Function to plot the trajectory of B relative to A
    plt.figure(figsize=(8, 6))
    plt.plot(y[:, 1], y[:, 0], label='Trajectory of B relative to A')
    plt.axis('equal')
    plt.xlabel('y (km)')
    plt.ylabel('x (km)')
    plt.title('Motion of chaser B relative to target A')
    plt.grid(True)
    plt.scatter(y[0, 1], y[0, 0], color='red', label="Start of B's trajectory")
    plt.legend()
    plt.show()

# Execute Example_7_03 function
Example_7_03()
