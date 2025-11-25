from LST import LST

# Example usage
if __name__ == "__main__":
    # Data declaration for Example 5.6:
    # East longitude:
    degrees = 139
    minutes = 47
    seconds = 0

    # Date:
    year = 2004
    month = 3
    day = 3

    # Universal time:
    hour = 4
    minute = 30
    second = 0

    # Convert negative (west) longitude to east longitude
    if degrees < 0:
        degrees += 360

    # Express the longitudes as decimal numbers
    EL = degrees + minutes / 60 + seconds / 3600
    WL = 360 - EL

    # Express universal time as a decimal number
    ut = hour + minute / 60 + second / 3600

    # Algorithm 5.3
    lst = LST(year, month, day, ut, EL)

    # Echo the input data and output the results to the console
    print("-----------------------------------------------------------")
    print("Example 5.6: Local sidereal time calculation")
    print("\nInput data:")
    print(f"Year = {year}")
    print(f"Month = {month}")
    print(f"Day = {day}")
    print(f"UT (hr) = {ut:.2f}")
    print(f"West Longitude (deg) = {WL:.2f}")
    print(f"East Longitude (deg) = {EL:.2f}\n")
    print("Solution:")
    print(f"Local Sidereal Time (deg) = {lst:.3f}")
    print(f"Local Sidereal Time (hr) = {lst / 15:.3f}")
    print("-----------------------------------------------------------")