import numpy as np
from ra_and_dec_from_r import ra_and_dec_from_r

r = np.array([-5368, -1784, 3691])
ra, dec = ra_and_dec_from_r(r)
print('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––')
print('\n Example 4.1')
print(f'\n r = [{r[0]} {r[1]} {r[2]}] (km)')
print(f'\n right ascension = {ra:.2f} deg')
print(f'\n declination = {dec:.2f} deg')
print('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')