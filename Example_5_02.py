import numpy as np
from coe_from_sv import coe_from_sv
from lambert import lambert

def example_5_02():
    deg = np.pi / 180
    mu = 398600  # Gravitational parameter (km^3/s^2)
    
    # Data declaration for Example 5.2:
    r1 = np.array([5000, 10000, 2100])
    r2 = np.array([-14600, 2500, 7000])
    dt = 3600  # Time between r1 and r2 (s)
    string = 'pro'
    
    # Solve Lambert's problem using the provided lambert function
    v1, v2 = lambert(r1, r2, dt, string)
    
    # Obtain orbital elements using r1 and v1
    coe = coe_from_sv(r1, v1, mu)
    
    # Save the initial true anomaly
    TA1 = coe[5]
    
    # Obtain orbital elements using r2 and v2
    coe = coe_from_sv(r2, v2, mu)
    
    # Save the final true anomaly
    TA2 = coe[5]
    
    # Echo the input data and output the results to the console
    print('-----------------------------------------------------------')
    print('Example 5.2: Lambert’s Problem')
    print('\nInput data:')
    print(f'Gravitational parameter (km^3/s^2) = {mu}')
    print(f'r1 (km) = [{r1[0]} {r1[1]} {r1[2]}]')
    print(f'r2 (km) = [{r2[0]} {r2[1]} {r2[2]}]')
    print(f'Elapsed time (s) = {dt}')
    
    print('\nSolution:')
    print(f'v1 (km/s) = [{v1[0]} {v1[1]} {v1[2]}]')
    print(f'v2 (km/s) = [{v2[0]} {v2[1]} {v2[2]}]')
    
    print('\nOrbital elements:')
    print(f'Angular momentum (km^2/s) = {coe[0]}')
    print(f'Eccentricity = {coe[1]}')
    print(f'Inclination (deg) = {coe[3] / deg}')
    print(f'RA of ascending node (deg) = {coe[2] / deg}')
    print(f'Argument of perigee (deg) = {coe[4] / deg}')
    print(f'True anomaly initial (deg) = {TA1 / deg}')
    print(f'True anomaly final (deg) = {TA2 / deg}')
    print(f'Semimajor axis (km) = {coe[6]}')
    print(f'Periapse radius (km) = {coe[0]**2 / mu / (1 + coe[1])}')
    
    # If the orbit is an ellipse, output its period
    if coe[1] < 1:
        T = 2 * np.pi / np.sqrt(mu) * coe[6]**1.5
        print('Period:')
        print(f'Seconds = {T}')
        print(f'Minutes = {T / 60}')
        print(f'Hours = {T / 3600}')
        print(f'Days = {T / 24 / 3600}')
    print('-----------------------------------------------------------')

# Run the example
if __name__ == "__main__":
    example_5_02()



