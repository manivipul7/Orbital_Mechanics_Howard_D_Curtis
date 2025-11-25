import numpy as np
from gibbs import gibbs
from coe_from_sv import coe_from_sv

def main():
    deg = np.pi / 180
    mu = 398600

    # Data declaration for Example 5.1
    r1 = np.array([-294.32, 4265.1, 5986.7])
    r2 = np.array([-1365.5, 3637.6, 6346.8])
    r3 = np.array([-2940.3, 2473.7, 6555.8])

    # Echo the input data
    print("--------------------------------------------------------")
    print("Example 5.1: Gibbs Method")
    print("\nInput data:")
    print(f"Gravitational parameter (km^3/s^2) = {mu}")
    print(f"r1 (km) = {r1}")
    print(f"r2 (km) = {r2}")
    print(f"r3 (km) = {r3}\n")

    # Algorithm 5.1
    v2, ierr = gibbs(r1, r2, r3, mu)

    # If the vectors r1, r2, r3 are not coplanar, abort
    if ierr == 1:
        print("These vectors are not coplanar.\n")
        return

    # Algorithm 4.2
    coe = coe_from_sv(r2, v2, mu)
    h, e, RA, incl, w, TA, a = coe

    # Output the results
    print("Solution:")
    print(f"v2 (km/s) = {v2}")
    print("\nOrbital elements:")
    print(f"Angular momentum (km^2/s) = {h}")
    print(f"Eccentricity = {e}")
    print(f"Inclination (deg) = {incl / deg}")
    print(f"RA of ascending node (deg) = {RA / deg}")
    print(f"Argument of perigee (deg) = {w / deg}")
    print(f"True anomaly (deg) = {TA / deg}")
    print(f"Semimajor axis (km) = {a}")

    # If the orbit is an ellipse, output the period
    if e < 1:
        T = 2 * np.pi / np.sqrt(mu) * a ** 1.5
        print(f"Period (s) = {T}")

    print("--------------------------------------------------------\n")

if __name__ == "__main__":
    main()