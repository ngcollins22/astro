# AOE 3154 
# Homework 2 - Nathan Collins

import numpy as np


print("Comet Orbit:")
#3
r = 6 #radius (AU)
v = 2 #velocity (AU/yr)
phi = -30 #flight path angle (deg)

#Calculate specific angular momentum
h = np.cos(phi * np.pi/180) * r * v #AU^2 /yr
print("#3: h =", f'{h:.3f}', "AU^2/yr")

mu = 39.48 #Gravitational parameter of sun (AU^3/yr^2)
eps = 0.5*np.square(v) - mu/r #Calculate specific energy
print("#4: specific energy =", f'{eps:.3f}')

e = np.sqrt(1 + 2*eps*np.square(h)/np.square(mu))

if eps < 0:
    if e == 0:
        print("#5: orbit is circular")
    else:
        print("#5: orbit is elliptic")
elif eps == 0:
    print("#5: orbit is parabolic")
else:
    print("#5: orbit is hyperbolic")

print("#6: eccentricity =", f'{e:.3f}')

p = np.square(h)/mu # Parameter of orbit
rp = p/(1+e) # perihelion radius of orbit

print("#7: perehelion radius =", f'{rp:.3f}')

a = rp/(1-e) #semimajor axis
TP = 2*np.pi*np.sqrt(a**3/mu) #orbital period

print("#8: orbital period =", f'{TP:.3f}')
print("-----")
print("Orbit of unkown object")

r = [1,2,4] # radius vector
v = [1,-4,-2] # velocity vector

h = np.cross(r, v) # specific angular momentum vector

print("#9: 1st component of specific angular momentum vector =", f'{h[0]:.3f}')
print("#10: 2nd =", f'{h[1]:.3f}')
print("#11: 3rd  =", f'{h[2]:.3f}')

e = np.cross(v,h)/mu - r/np.linalg.norm(r)

print("#12: 1st component of eccentricity vector =", f'{e[0]:.3f}')
print("#13: 2nd =", f'{e[1]:.3f}')
print("#14: 3rd  =", f'{e[2]:.3f}')

e_scalar = np.linalg.norm(e)

if e_scalar == 0:
    print("#15: circular, less than 0, enter circular trajectory")
elif e_scalar < 1:
    print("#15: elliptical, less than 0, return again")
elif e_scalar == 1:
    print("#15: parabolic, equal to 0, fly off into space")
else:
    print("#15: hyperbolic, greater than 0, fly off into space")

mu = 39.48 #grav param of sun
v_scalar = np.linalg.norm(v) # calc magnitudes
r_scalar = np.linalg.norm(r)

v_inf = np.sqrt(v_scalar**2 - 2*mu/r_scalar) #calc excess energy

print("#16: excess velocity =", f'{v_inf:.3f}')

delta = 2*np.arcsin(1/e_scalar) * 180/np.pi # calc turning angle

print("#17: turning angle =", f'{delta:.3f}')



