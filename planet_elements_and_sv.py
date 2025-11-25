import numpy as np
from math import radians, degrees, sqrt, atan, sin, cos
from kepler_E import kepler_E
from J0 import J0
from sv_from_coe import sv_from_coe

# Define global constants
mu = 1.32712440018e11  # gravitational parameter of the sun (km^3/s^2)
deg = np.pi / 180  # conversion factor between degrees and radians

def planet_elements_and_sv(planet_id, year, month, day, hour, minute, second):
    # Function to calculate orbital elements and state vector of a planet
    
    # Convert date and time to Julian day
    j0 = J0(year, month, day)
    ut = (hour + minute / 60 + second / 3600) / 24
    jd = j0 + ut
    
    # Load J2000 orbital elements and rates from Table 8.1
    J2000_coe, rates = planetary_elements(planet_id)
    
    # Calculate centuries between J2000 and the given epoch
    t0 = (jd - 2451545) / 36525
    
    # Calculate elements at the given epoch
    elements = J2000_coe + rates * t0
    
    # Extract elements from the array
    a = elements[0]
    e = elements[1]

    h = sqrt(mu * a * (1 - e**2))

    incl = elements[2]
    RA = zero_to_360(elements[3])
    w_hat = zero_to_360(elements[4])
    L = zero_to_360(elements[5])
    
    # Calculate derived elements
    w = zero_to_360(w_hat - RA)
    M = zero_to_360(L - w_hat)
    
    # Calculate eccentric anomaly E
    E = kepler_E(e, M*deg)    
    
    # Calculate true anomaly TA
    TA = zero_to_360(2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2)) / deg)
    
    # Formulate orbital elements vector (converted radians to degrees for RA, incl, w)
    coe = [h, e, RA, incl, w, TA, a, w_hat, L, M, E/deg]
    
    # Calculate state vector (position and velocity) using the orbital elements
    r, v = sv_from_coe([h, e, RA*deg, incl*deg, w*deg, TA*deg], mu)
    
    return coe, r, v, jd

def planetary_elements(planet_id):
    # Function to extract J2000 orbital elements and centennial rates for a planet
    
    # J2000 elements and rates from Table 8.1 (converted to km and degrees)
    J2000_elements = np.array([
        [0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084],
        [0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973],
        [1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 100.46435],
        [1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332],
        [5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438],
        [9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432],
        [19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218],
        [30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003],
        [39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881]
    ])
    
    cent_rates = np.array([
      [ 0.00000066, 0.00002527, -23.51, -446.30, 573.57, 538101628.29],
        [0.00000092, -0.00004938, -2.86, -996.89, -108.80, 210664136.06],
        [-0.00000005, -0.00003804, -46.94, -18228.25, 1198.28, 129597740.63],
        [-0.00007221, 0.00011902, -25.47, -1020.19, 1560.78, 68905103.78],
        [0.00060737, -0.00012880, -4.15, 1217.17, 839.93, 10925078.35],
        [-0.00301530, -0.00036762, 6.11, -1591.05, -1948.89,  4401052.95],
        [0.00152025, -0.00019150, -2.09, -1681.4, 1312.56, 1542547.79],
        [-0.00125196, 0.00002514, -3.64, -151.25, -844.43, 786449.21],
        [-0.00076912, 0.00006465, 11.07, -37.33, -132.25, 522747.90]
    ])
    
    # Extract J2000 elements and rates for the given planet_id
    J2000_coe = J2000_elements[planet_id - 1,:]
    rates = cent_rates[planet_id - 1,:]

    # Convert AU to km and arcseconds to degrees
    au = 149597871  # astronomical unit (km)
    J2000_coe[0] *= au
    rates[0] *= au

    # Convert arcseconds to degrees
    rates[2:6] /= 3600
    
    return J2000_coe, rates

def zero_to_360(x):
    # Function to reduce an angle to lie within 0-360 degrees
    if x >= 360:
        x = x - np.floor(x / 360) * 360
    elif x < 0:
        x = x - (np.floor(x / 360) - 1) * 360
    return x % 360
