import numpy as np
import math
from coe_from_sv import coe_from_sv
from interplanetary import interplanetary
from month_planet_names import month_planet_names

# Global constants
mu = 1.327124e11  # gravitational parameter of the sun (km^3/s^2)
deg = math.pi / 180

# Data declaration for Example 8.8:
# Departure
depart_planet_id = 3
depart_year = 1996
depart_month = 11
depart_day = 7
depart_hour = 0
depart_minute = 0
depart_second = 0
depart = [depart_planet_id, depart_year, depart_month, depart_day, depart_hour, depart_minute, depart_second]

# Arrival
arrive_planet_id = 4
arrive_year = 1997
arrive_month = 9
arrive_day = 12
arrive_hour = 0
arrive_minute = 0
arrive_second = 0
arrive = [arrive_planet_id, arrive_year, arrive_month, arrive_day, arrive_hour, arrive_minute, arrive_second]

# Algorithm 8.2:
planet1, planet2, trajectory = interplanetary(depart, arrive)
print(planet1)
print(planet2)

# Assuming planet1[0] is a list or array, correct the slicing method
R1 = np.array(planet1[0])
Vp1 = np.array(planet1[1])
jd1 = planet1[2]

# Similarly for planet2
R2 = np.array(planet2[0])
Vp2 = np.array(planet2[1])
jd2 = planet2[2]

V1 = np.array(trajectory[0])
V2 = np.array(trajectory[1])

tof = jd2 - jd1

# Use Algorithm 4.2 to find the orbital elements of the spacecraft trajectory based on [Rp1, V1]...
coe = coe_from_sv(R1, V1, mu)
# ... and [R2, V2]
coe2 = coe_from_sv(R2, V2, mu)

# Equations 8.94 and 8.95:
vinf1 = V1 - Vp1
vinf2 = V2 - Vp2

# Echo the input data and output the solution to the command window:
def print_results():
    print("------------------------------------------------------------")
    print("Example 8.8")
    print("\nDeparture:\n")
    month_name_depart, planet_name_depart = month_planet_names(depart[2], depart[0])
    print(f"Planet: {planet_name_depart}")
    print(f"Year: {depart[1]}")
    print(f"Month: {month_name_depart}")
    print(f"Day: {depart[3]}")
    print(f"Hour: {depart[4]}")
    print(f"Minute: {depart[5]}")
    print(f"Second: {depart[6]}")
    print(f"\nJulian day: {jd1:.3f}\n")
    print(f"Planet position vector (km): [{R1[0]} {R1[1]} {R1[2]}]")
    print(f"Magnitude: {np.linalg.norm(R1):.3f}")
    print(f"Planet velocity (km/s): [{Vp1[0]} {Vp1[1]} {Vp1[2]}]")
    print(f"Magnitude: {np.linalg.norm(Vp1):.3f}")
    print(f"Spacecraft velocity (km/s): [{V1[0]} {V1[1]} {V1[2]}]")
    print(f"Magnitude: {np.linalg.norm(V1):.3f}")
    print(f"v-infinity at departure (km/s): [{vinf1[0]} {vinf1[1]} {vinf1[2]}]")
    print(f"Magnitude: {np.linalg.norm(vinf1):.3f}")
    print(f"\nTime of flight: {tof} days\n")
    print("\nArrival:\n")
    month_name_arrive, planet_name_arrive = month_planet_names(arrive[2], arrive[0])
    print(f"Planet: {planet_name_arrive}")
    print(f"Year: {arrive[1]}")
    print(f"Month: {month_name_arrive}")
    print(f"Day: {arrive[3]}")
    print(f"Hour: {arrive[4]}")
    print(f"Minute: {arrive[5]}")
    print(f"Second: {arrive[6]}")
    print(f"\nJulian day: {jd2:.3f}\n")
    print(f"Planet position vector (km): [{R2[0]} {R2[1]} {R2[2]}]")
    print(f"Magnitude: {np.linalg.norm(R2):.3f}")
    print(f"Planet velocity (km/s): [{Vp2[0]} {Vp2[1]} {Vp2[2]}]")
    print(f"Magnitude: {np.linalg.norm(Vp2):.3f}")
    print(f"Spacecraft velocity (km/s): [{V2[0]} {V2[1]} {V2[2]}]")
    print(f"Magnitude: {np.linalg.norm(V2):.3f}")
    print(f"v-infinity at arrival (km/s): [{vinf2[0]} {vinf2[1]} {vinf2[2]}]")
    print(f"Magnitude: {np.linalg.norm(vinf2):.3f}")
    print("\n\nOrbital elements of flight trajectory:\n")
    print(f"Angular momentum (km^2/s): {coe[0]}")
    print(f"Eccentricity: {coe[1]}")
    print(f"Right ascension of the ascending node (deg): {coe[2] / deg}")
    print(f"Inclination to the ecliptic (deg): {coe[3] / deg}")
    print(f"Argument of perihelion (deg): {coe[4] / deg}")
    print(f"True anomaly at departure (deg): {coe[5] / deg}")
    print(f"True anomaly at arrival (deg): {coe2[5] / deg}")
    print(f"Semimajor axis (km): {coe[6]}")
    if coe[1] < 1:  # If the orbit is an ellipse, output the period
        period = 2 * math.pi / math.sqrt(mu) * coe[6] ** 1.5 / (24 * 3600)
        print(f"Period (days): {period}")
    print("------------------------------------------------------------")

print_results()
