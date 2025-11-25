""" A system consists of two massive bodies m1 and m2 each having a mass of 10^26 kg. At time t ¼ 0 the state vectors of the two
particles in an inertial frame are r1 = 0 , r2 = 3000 i km, v1 = 10 i + 20 j + 30 k km/s, and v2 = 40 j km/s.
The function twobody3d.py in Appendix D.5 contains within it the data for this problem. Embedded in the
program is the subfunction rates, which computes the accelerations given by Eqs. (2.19a) and (2.19b). twobody3d.m
uses the solution vector from rkf45.m to plot Figs. 2.2 and 2.3, which summarize the results requested in the problem
statement.
In answer to part (a), Fig. 2.2 shows the motion of the two-body system relative to the inertial frame. m1 and m2 are soon
established in a periodic helical motion around the straight-line trajectory of the center of mass G through space. This
pattern continues indefinitely.
Fig. 2.3(a) relates to part (b) of the problem. The very same motion appears rather less complex when viewed from m1.
In fact, we see that R2(t)  R1(t), the trajectory of m2 relative to m1, appears to be an elliptical path. So does RG(t)  R1(t),
the path of the center of mass around m1.
Finally, for part (c) of the problem, Fig. 2.3(b) reveals that both m1 and m2 follow apparently elliptical paths around the
center of mass
"""

# Define the gravitational constant
G = 6.6742e-20  # km^3/kg/s^2

# Define the masses of the two bodies
m1 = 1e26  # kg
m2 = 1e26
M = m1 + m2  # kg

# Define the initial positions and velocities of the two bodies
r1_0 = np.array([0, 0, 0])  # km
v1_0 = np.array([10, 20, 30])  # km/s
r2_0 = np.array([3000, 0, 0])  # km
v2_0 = np.array([0, 40, 0])  # km/s

# Define the initial state vector
f0 = np.array([r1_0, v1_0, r2_0, v2_0])

# Define the rates function




