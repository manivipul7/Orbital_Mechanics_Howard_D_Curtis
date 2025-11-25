import numpy as np
import matplotlib.pyplot as plt
from heun import heun

# Define the differential equation
def rates(t, f):
    x = f[0]
    Dx = f[1]
    D2x = (Fo / m) * np.sin(w * t) - 2 * z * wn * Dx - wn**2 * x
    return np.array([Dx, D2x])

# System properties
m = 1
z = 0.03
wn = 1
Fo = 1
w = 0.4 * wn

# Time range
t0 = 0
tf = 110
tspan = [t0, tf]

# Initial conditions
x0 = 0
Dx0 = 0
f0 = np.array([x0, Dx0])

# Calculate and plot the solution for h = 1.0
h1 = 1.0
t1, f1 = heun(rates, tspan, f0, h1)

# Calculate and plot the solution for h = 0.1
h2 = 0.1
t2, f2 = heun(rates, tspan, f0, h2)

# Plot solutions
plt.figure(figsize=(10, 6))

plt.plot(t1, f1[:, 0], '-r', linewidth=0.5, label='h = 1.0')
plt.plot(t2, f2[:, 0], '-k', linewidth=1, label='h = 0.1')

plt.xlabel('time, s')
plt.ylabel('x, m')
plt.grid(True)
plt.axis([0, 110, -2, 2])
plt.legend()
plt.title('Response of a Damped Single Degree of Freedom Spring-Mass System')
plt.show()
