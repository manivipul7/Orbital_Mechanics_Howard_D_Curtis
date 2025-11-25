import numpy as np
import matplotlib.pyplot as plt
from rkf1_4 import rkf1_4

# Constants and parameters
m = 1
z = 0.03
wn = 1
Fo = 1
w = 0.4 * wn
x0 = 0
x_dot0 = 0
f0 = np.array([x0, x_dot0])
t0 = 0
tf = 110
tspan = [t0, tf]

# Define the rates function
def rates(t, f):
    x = f[0]
    Dx = f[1]
    D2x = (Fo / m) * np.sin(w * t) - 2 * z * wn * Dx - wn**2 * x
    return np.array([Dx, D2x])

# Solve using RK1 through RK4 with different time steps
results = {}
for rk in range(1, 5):
    if rk == 1:
        hs = [0.01, 0.1]
    elif rk == 2:
        hs = [0.1, 0.5]
    elif rk == 3:
        hs = [0.5, 1.0]
    elif rk == 4:
        hs = [1.0, 2.0]
    
    for h in hs:
        t, f = rkf1_4(rates, tspan, f0, h, rk)
        results[(rk, h)] = (t, f)

# Exact solution
wd = wn * np.sqrt(1 - z**2)
den = (wn**2 - w**2)**2 + (2 * w * wn * z)**2
C1 = (wn**2 - w**2) / den * Fo / m
C2 = -2 * w * wn * z / den * Fo / m
A = (x0 * wn / wd + x_dot0 / wd +
     (w**2 + (2 * z**2 - 1) * wn**2) / den * w / wd * Fo / m)
B = x0 + 2 * w * wn * z / den * Fo / m
t_exact = np.linspace(t0, tf, 5000)
x_exact = (A * np.sin(wd * t_exact) + B * np.cos(wd * t_exact)) * np.exp(-wn * z * t_exact) + C1 * np.sin(w * t_exact) + C2 * np.cos(w * t_exact)

# Plot solutions
plt.figure(figsize=(10, 15))

# Plot exact solution
plt.subplot(5, 1, 1)
plt.plot(t_exact / max(t_exact), x_exact / max(x_exact), 'k', linewidth=1)
plt.grid(False)
plt.axis('tight')
plt.title('Exact')

# Plot RK1 solutions
plt.subplot(5, 1, 2)
for h in [0.01, 0.1]:
    t, f = results[(1, h)]
    plt.plot(t / max(t), f[:, 0] / max(f[:, 0]), label=f'h = {h}', linewidth=1)
plt.grid(False)
plt.axis('tight')
plt.title('RK1')
plt.legend()

# Plot RK2 solutions
plt.subplot(5, 1, 3)
for h in [0.1, 0.5]:
    t, f = results[(2, h)]
    plt.plot(t / max(t), f[:, 0] / max(f[:, 0]), label=f'h = {h}', linewidth=1)
plt.grid(False)
plt.axis('tight')
plt.title('RK2')
plt.legend()

# Plot RK3 solutions
plt.subplot(5, 1, 4)
for h in [0.5, 1.0]:
    t, f = results[(3, h)]
    plt.plot(t / max(t), f[:, 0] / max(f[:, 0]), label=f'h = {h}', linewidth=1)
plt.grid(False)
plt.axis('tight')
plt.title('RK3')
plt.legend()

# Plot RK4 solutions
plt.subplot(5, 1, 5)
for h in [1.0, 2.0]:
    t, f = results[(4, h)]
    plt.plot(t / max(t), f[:, 0] / max(f[:, 0]), label=f'h = {h}', linewidth=1)
plt.grid(False)
plt.axis('tight')
plt.title('RK4')
plt.legend()

plt.tight_layout()
plt.show()
