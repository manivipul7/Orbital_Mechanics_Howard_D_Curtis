from kepler_E import kepler_E

# Example data from Example 3.2
e = 0.37255
M = 3.6029

# Compute eccentric anomaly E using kepler_E function
E = kepler_E(e, M)

# Print results to the console
print('–––––––––––––––––––––––––––––––––––––––––––––––––––––––')
print('Example 3.2')
print(f'\nEccentricity: {e}')
print(f'Mean anomaly (radians): {M}')
print(f'Eccentric anomaly (radians): {E}')
print('–––––––––––––––––––––––––––––––––––––––––––––––––––––––')