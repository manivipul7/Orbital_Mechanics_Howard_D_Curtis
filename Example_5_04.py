from J0 import J0

def example_5_04():
    """
    This program computes J0 and the Julian day number using the data
    in Example 5.4.
    
    Parameters:
    None
    
    Returns:
    None
    """
    # Data declaration for Example 5.4:
    year = 2004
    month = 5
    day = 12
    hour = 14
    minute = 45
    second = 30
    
    # Calculate universal time in hours
    ut = hour + minute / 60 + second / 3600
    
    # Equation 5.46
    j0 = J0(year, month, day)
    
    # Equation 5.47
    jd = j0 + ut / 24
    
    # Echo the input data and output the results to the console
    print('-----------------------------------------------------------')
    print('Example 5.4: Julian day calculation')
    print('\nInput data:')
    print(f'\nYear = {year}')
    print(f'\nMonth = {month}')
    print(f'\nDay = {day}')
    print(f'\nHour = {hour}')
    print(f'\nMinute = {minute}')
    print(f'\nSecond = {second}')
    print(f'\nJulian day number = {jd:.3f}')
    print('-----------------------------------------------------------')

# Run the example
if __name__ == "__main__":
    example_5_04()