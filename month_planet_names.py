def month_planet_names(month_id, planet_id):
    # Define lists containing names of months and planets
    months = ['January', 'February', 'March', 'April', 'May', 'June',
              'July', 'August', 'September', 'October', 'November', 'December']
    
    planets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter',
               'Saturn', 'Uranus', 'Neptune', 'Pluto']
    
    # Adjust month_id and planet_id to be zero-indexed for list access
    month_index = month_id - 1
    planet_index = planet_id - 1
    
    # Check if the given IDs are within valid range
    if 0 <= month_index < len(months):
        month_name = months[month_index]
    else:
        month_name = 'Invalid month ID'
    
    if 0 <= planet_index < len(planets):
        planet_name = planets[planet_index]
    else:
        planet_name = 'Invalid planet ID'
    
    return month_name, planet_name
