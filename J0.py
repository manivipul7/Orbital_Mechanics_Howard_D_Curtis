def J0(year, month, day):
    """
    This function computes the Julian day number at 0 UT for any year
    between 1900 and 2100 using Equation 5.48.
    
    Parameters:
    year (int): Year (range: 1901-2099)
    month (int): Month (range: 1-12)
    day (int): Day (range: 1-31)
    
    Returns:
    float: Julian day at 0 hr UT (Universal Time)
    """
    j0 = 367 * year - (7 * (year + ((month + 9) // 12)) // 4) + (275 * month // 9) + day + 1721013.5
    return j0
