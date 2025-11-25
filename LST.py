from J0 import J0

def zero_to_360(x):
    """
    This function reduces an angle to the range 0-360 degrees.

    Parameters:
    x (float): The angle (degrees) to be reduced

    Returns:
    float: The reduced value
    """
    if x >= 360:
        x = x - (x // 360) * 360
    elif x < 0:
        x = x - (int(x / 360) - 1) * 360
    return x

def LST(y, m, d, ut, EL):
    """
    This function calculates the local sidereal time.

    Parameters:
    y (int): Year
    m (int): Month
    d (int): Day
    ut (float): Universal Time (hours)
    EL (float): East longitude (degrees)

    Returns:
    float: Local sidereal time (degrees)
    """
    # Equation 5.48
    j0 = J0(y, m, d)

    # Equation 5.49
    j = (j0 - 2451545) / 36525

    # Equation 5.50
    g0 = 100.4606184 + 36000.77004 * j + 0.000387933 * j**2 - 2.583e-8 * j**3

    # Reduce g0 to the range 0-360 degrees
    g0 = zero_to_360(g0)

    # Equation 5.51
    gst = g0 + 360.98564724 * ut / 24

    # Equation 5.52
    lst = gst + EL

    # Reduce lst to the range 0-360 degrees
    lst = zero_to_360(lst)

    return lst


