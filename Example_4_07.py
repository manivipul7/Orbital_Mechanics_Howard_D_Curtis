import numpy as np
from sv_from_coe import sv_from_coe

# Example usage
def example_4_07():
    deg = np.pi / 180
    mu = 398600

    # Data declaration for Example 4.7 (angles in degrees)
    h = 80000
    e = 1.4
    RA = 40
    incl = 30
    w = 60
    TA = 30

    coe = [h, e, RA * deg, incl * deg, w * deg, TA * deg]

    # Algorithm 4.5 (requires angular elements be in radians)
    r, v = sv_from_coe(coe, mu)

    # Output the results
    print("–––––––––––––––––––––––––––––––––––––––––––––––––––––")
    print("Example 4.7")
    print(f"\nGravitational parameter (km^3/s^2) = {mu}")
    print(f"\nAngular momentum (km^2/s) = {h}")
    print(f"\nEccentricity = {e}")
    print(f"\nRight ascension (deg) = {RA}")
    print(f"\nArgument of perigee (deg) = {w}")
    print(f"\nTrue anomaly (deg) = {TA}")
    print(f"\nState vector:")
    print(f"r (km) = [{r[0]}, {r[1]}, {r[2]}]")
    print(f"v (km/s) = [{v[0]}, {v[1]}, {v[2]}]")
    print("–––––––––––––––––––––––––––––––––––––––––––––––––––––")

example_4_07()
