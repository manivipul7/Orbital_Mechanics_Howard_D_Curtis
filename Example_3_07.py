import numpy as np
from rv_from_r0v0 import rv_from_r0v0

# Data declaration for Example 3.7
R0 = np.array([7000, -12124, 0])  # initial position vector (example values in km)
V0 = np.array([2.6679, 4.6210, 0])  # initial velocity vector (example values in km/s)
t = 3600  # elapsed time (example value in seconds)
# Algorithm 3.4
R, V = rv_from_r0v0(R0, V0, t)
# Echo the input data and output the results
print('–––––––––––––––––––––––––––––––––––––––––––––––––––––')
print(' Example 3.7')
print('\n Initial position vector (km):')
print(f' r0 = ({R0[0]}, {R0[1]}, {R0[2]})')
print('\n Initial velocity vector (km/s):')
print(f' v0 = ({V0[0]}, {V0[1]}, {V0[2]})')
print(f'\n\n Elapsed time = {t} s')
print('\n Final position vector (km):')
print(f' r = ({R[0]}, {R[1]}, {R[2]})')
print('\n Final velocity vector (km/s):')
print(f' v = ({V[0]}, {V[1]}, {V[2]})')
print('–––––––––––––––––––––––––––––––––––––––––––––––––––––')