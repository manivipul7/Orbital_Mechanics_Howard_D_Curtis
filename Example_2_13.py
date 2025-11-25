import numpy as np
from rv_from_r0v0_ta import rv_from_r0v0_ta

def example_2_13():
    mu = 398600
    
    # Input data
    R0 = np.array([8182.4, -6865.9, 0])  # Initial position vector (km)
    V0 = np.array([0.47572, 8.8116, 0])  # Initial velocity vector (km/s)
    dt = 120  # Change in true anomaly (degrees)
    
    # Algorithm 2.3
    R, V = rv_from_r0v0_ta(R0, V0, dt, mu)
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)
    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)
    
    # Output results
    print('––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––')
    print('Example 2.13')
    print('\nInitial state vector:')
    print(f'\nr = [{R0[0]:.4f}, {R0[1]:.4f}, {R0[2]:.4f}] (km)')
    print(f'Magnitude = {r0:.4f} km')
    print(f'\nv = [{V0[0]:.4f}, {V0[1]:.4f}, {V0[2]:.4f}] (km/s)')
    print(f'Magnitude = {v0:.4f} km/s')
    print(f'\nState vector after {dt} degree change in true anomaly:')
    print(f'\nr = [{R[0]:.4f}, {R[1]:.4f}, {R[2]:.4f}] (km)')
    print(f'Magnitude = {r:.4f} km')
    print(f'\nv = [{V[0]:.4f}, {V[1]:.4f}, {V[2]:.4f}] (km/s)')
    print(f'Magnitude = {v:.4f} km/s')
    print('––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––')

if __name__ == "__main__":
    example_2_13()
