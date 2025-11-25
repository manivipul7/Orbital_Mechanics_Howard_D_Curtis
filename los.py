import numpy as np

def los(r_sat, r_sun):
    """
    This function uses the ECI position vectors of the satellite (r_sat)
    and the sun (r_sun) to determine whether the Earth is in the line of
    sight between the two.

    Parameters:
        r_sat (numpy array): ECI position vector of the satellite
        r_sun (numpy array): ECI position vector of the sun

    Returns:
        light_switch (int): 0 if the Earth is in the line of sight, 1 otherwise
    """
    # Earth's radius (km)
    RE = 6378

    # Magnitudes of the position vectors
    rsat = np.linalg.norm(r_sat)
    rsun = np.linalg.norm(r_sun)

    # Angle between sun and satellite position vectors
    theta = np.degrees(np.arccos(np.dot(r_sat, r_sun) / (rsat * rsun)))

    # Angle between the satellite position vector and the radial to the point
    # of tangency with the Earth of a line from the satellite
    theta_sat = np.degrees(np.arccos(RE / rsat))

    # Angle between the sun position vector and the radial to the point
    # of tangency with the Earth of a line from the sun
    theta_sun = np.degrees(np.arccos(RE / rsun))

    # Determine whether a line from the sun to the satellite intersects the Earth
    if theta_sat + theta_sun <= theta:
        light_switch = 0  # Yes
    else:
        light_switch = 1  # No

    return light_switch

