# Class from 9/11/24

import numpy as np

# Generic planet
Rp = 1 # DU (radius of planet)

#add massless mountain
alt = 0.1 #DU
r = Rp + alt # DU to top of mountain

v = 4 # DU/TU

# does the cannonball hit planet?
# if radius of periapsis does, then it'll hit the planet
# use energy equation

mu = 1 # in canonical units
E = (v**2)/2 - mu/r
print(E)
# hyperbolic orbit, thing is defintely not hitting the planet

v = 0.1
E = (v**2)/2 - mu/r
print(E)
# negative specific energy, closed orbit
a = -1



