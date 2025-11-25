import numpy as np
import matplotlib.pyplot as plt
from sv_from_coe import sv_from_coe
from rva_relative import rva_relative
from rv_from_r0v0 import rv_from_r0v0

# Constants
mu = 398600  # Gravitational parameter (km^3/s^2)
RE = 6378  # Earth's radius (km)
deg = np.pi / 180  # Degree to radian conversion

# Orbital parameters for Spacecraft A and B
h_A = 52059
e_A = 0.025724
i_A = 60 * deg
RAAN_A = 40 * deg
omega_A = 30 * deg
theta_A = 40 * deg

h_B = 52362
e_B = 0.0072696
i_B = 50 * deg
RAAN_B = 40 * deg
omega_B = 120 * deg
theta_B = 40 * deg

vdir = np.array([1, 1, 1])

# Compute state vectors for Spacecraft A and B
coe_A = [h_A, e_A, RAAN_A, i_A, omega_A, theta_A]
coe_B = [h_B, e_B, RAAN_B, i_B, omega_B, theta_B]

rA0, vA0 = sv_from_coe(coe_A, mu)
rB0, vB0 = sv_from_coe(coe_B, mu)

h0 = np.cross(rA0, vA0)

# Period of spacecraft A
TA = 2 * np.pi / mu**2 * (h_A/np.sqrt(1 - e_A**2))**3

# Number of time steps per period of A's orbit
n = 100

# Time step as a fraction of A's period
dt = TA / n

# Number of periods of A's orbit for which the trajectory will be plotted
n_Periods = 60

# Initialize time
t = -dt

x = np.zeros(n_Periods * n)
y = np.zeros(n_Periods * n)
z = np.zeros(n_Periods * n)
r = np.zeros(n_Periods * n)
T = np.zeros(n_Periods * n)

# Generate the trajectory of B relative to A
for count in range(n_Periods * n):
    # Update time
    t += dt
    
    # Update state vectors of A and B
    rA, vA = rv_from_r0v0(rA0, vA0, t)
    rB, vB = rv_from_r0v0(rB0, vB0, t)
    
    # Compute relative position using Algorithm 7.1
    r_rel, v_rel, a_rel = rva_relative(rA, vA, rB, vB, mu)
    
    # Store components of relative position vector
    x[count] = r_rel[0]
    y[count] = r_rel[1]
    z[count] = r_rel[2]
    r[count] = np.linalg.norm(r_rel)
    T[count] = t

# Plot the trajectory of B relative to A
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z)
ax.axis('equal')
ax.grid(True)
# Make the panes (which form the sides of the 3D box) transparent
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
# Optionally, also make the grid lines transparent
ax.xaxis.pane.set_alpha(0)
ax.yaxis.pane.set_alpha(0)
ax.zaxis.pane.set_alpha(0)

vdir = [30,45]
ax.view_init(elev=vdir[0], azim=vdir[1])

# Draw co-moving x, y, z axes
ax.plot([0, 4000], [0, 0], [0, 0], color='black')  # x-axis
ax.text(4000, 0, 0, 'x', color='black')
ax.plot([0, 0], [0, 7000], [0, 0], color='black')  # y-axis
ax.text(0, 7000, 0, 'y', color='black')
ax.plot([0, 0], [0, 0], [0, 4000], color='black')  # z-axis
ax.text(0, 0, 4000, 'z', color='black')

# Label the origin of the moving frame attached to A
ax.text(0, 0, 0, 'A', color='black')

# Label the start of B's relative trajectory
ax.text(x[0], y[0], z[0], 'B', color='black')

# Draw the initial position vector of B
ax.plot([0, x[0]], [0, y[0]], [0, z[0]], color='red')

# Show plot
plt.show()