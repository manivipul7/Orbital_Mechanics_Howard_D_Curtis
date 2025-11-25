from lambert import lambert
from planet_elements_and_sv import planet_elements_and_sv
from izzo2015 import izzo2015
import numpy as np

def interplanetary(depart, arrive):
    """
    This function determines the spacecraft trajectory from the sphere
    of influence of planet 1 to that of planet 2 using Algorithm 8.2.
    """
    global mu

    # Extract departure data
    planet_id = depart[0]
    year = depart[1]
    month = depart[2]
    day = depart[3]
    hour = depart[4]
    minute = depart[5]
    second = depart[6]

    # Use Algorithm 8.1 to obtain planet 1’s state vector (don’t need its orbital elements)
    dum, Rp1, Vp1, jd1 = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second)

    # Extract arrival data
    planet_id = arrive[0]
    year = arrive[1]
    month = arrive[2]
    day = arrive[3]
    hour = arrive[4]
    minute = arrive[5]
    second = arrive[6]

    # Use Algorithm 8.1 to obtain planet 2’s state vector
    dum, Rp2, Vp2, jd2 = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second)

    tof = (jd2 - jd1) * 24 * 3600  # Time of flight from planet 1 to planet 2 in seconds

    # Patched conic assumption
    R1 = Rp1
    R2 = Rp2

    # Use Algorithm 5.2 to find the spacecraft’s velocity at departure and arrival, assuming a prograde trajectory
    V1, V2 = izzo2015(mu, R1, R2, tof, prograde=True)

    planet1 = [Rp1, Vp1, jd1]
    planet2 = [Rp2, Vp2, jd2]
    trajectory = [V1, V2]

    return planet1, planet2, trajectory

# Global constant
mu = 1.327124e11  # gravitational parameter of the sun (km^3/s^2)