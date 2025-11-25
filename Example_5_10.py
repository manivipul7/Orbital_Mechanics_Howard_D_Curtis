from rv_from_observe import rv_from_observe
from coe_from_sv import coe_from_sv
import numpy as np

def example_5_10():
    # Constants
    deg = np.pi / 180
    f = 1 / 298.256421867
    Re = 6378.13655
    wE = 7.292115e-5
    mu = 398600.4418

    # Data declaration for Example 5.10
    rho = 2551
    rhodot = 0
    A = 90
    Adot = 0.1130
    a = 30
    adot = 0.05651
    theta = 300
    phi = 60
    H = 0

    # Algorithm 5.4
    r, v = rv_from_observe(rho, rhodot, A, Adot, a, adot, theta, phi, H, Re, f, wE)

    # Algorithm 4.2
    coe = coe_from_sv(r, v, mu)
    h = coe[0]
    e = coe[1]
    RA = coe[2]
    incl = coe[3]
    w = coe[4]
    TA = coe[5]
    a = coe[6]

    # Equation 2.40
    rp = h**2 / mu / (1 + e)

    # Output the solution
    print('–––––––––––––––––––––––––––––––––––––––––––––––––––––')
    print('Example 5.10')
    print('\nInput data:\n')
    print(f'Slant range (km) = {rho}')
    print(f'Slant range rate (km/s) = {rhodot}')
    print(f'Azimuth (deg) = {A}')
    print(f'Azimuth rate (deg/s) = {Adot}')
    print(f'Elevation (deg) = {a}')
    print(f'Elevation rate (deg/s) = {adot}')
    print(f'Local sidereal time (deg) = {theta}')
    print(f'Latitude (deg) = {phi}')
    print(f'Altitude above sea level (km) = {H}')
    print('\nSolution:')
    print('\nState vector:\n')
    print(f'r (km) = [{r[0]}, {r[1]}, {r[2]}]')
    print(f'v (km/s) = [{v[0]}, {v[1]}, {v[2]}]')
    print('\nOrbital elements:\n')
    print(f'Angular momentum (km^2/s) = {h}')
    print(f'Eccentricity = {e}')
    print(f'Inclination (deg) = {incl / deg}')
    print(f'RA of ascending node (deg) = {RA / deg}')
    print(f'Argument of perigee (deg) = {w / deg}')
    print(f'True anomaly (deg) = {TA / deg}')
    print(f'Semimajor axis (km) = {a}')
    print(f'Perigee radius (km) = {rp}')
    if e < 1:
        T = 2 * np.pi / np.sqrt(mu) * a**1.5
        print(f'\nPeriod:')
        print(f'Seconds = {T}')
        print(f'Minutes = {T / 60}')
        print(f'Hours = {T / 3600}')
        print(f'Days = {T / 24 / 3600}')
    print('–––––––––––––––––––––––––––––––––––––––––––––––––––––')

# Run the example
example_5_10()