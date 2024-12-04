# AOE 3154 
# Homework 10 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u


#2

rpark = 6700 #km
rfinal = 7100 #km
# Constellation is 97deg:24/4/1
#find delta raan between planes in this constellation
[deltaRAAN, _, _] = u.walker(24,4,2,1)

#find mean motion of parking and final orbit
nbarpark = u.mean_motion(rpark, u.MU_EARTH)
nbarfinal = u.mean_motion(rfinal, u.MU_EARTH)

#find regression of nodes for each
raandot_final = u.RAAN_dot(rfinal, nbarfinal, 97*np.pi/180)
raandot_park = u.RAAN_dot(rpark, nbarpark, 97*np.pi/180)

#find relative regression of nodes
deltaraandot = raandot_park - raandot_final
#print(deltaRAAN)

#divide to find time to relative-regress the right amount
deltaT = deltaRAAN / deltaraandot

print("#2 Time to wait (days): ", f'{deltaT*(1/3600)*(1/24):.4f}')

#4/5
mu_sun = 1.327*10**11
mu_moon = 4.905*10**3
r_earth = 1 #AU

#you can use mu or mass of bodies - units cancel completely as well
rsoi_earth = u.rsoi(r_earth, u.MU_EARTH, mu_sun)
print("#4 Radius of sphere of influence of earth (AU): ", f'{rsoi_earth:.4f}')

frac_soi_moon = (mu_moon/u.MU_EARTH)**(2/5)
print("#5 Radius of sphere of influence of moon (fraction of earth-moon distance): ", f'{frac_soi_moon:.4f}')


mu_sun = 39.4784 #Canonical units
re = 1
rm = 1.524 # AU
at = (re + rm)/2

TOF = u.timePeriod(at, mu_sun)/2

et = (rm - re)/(rm + re)
pt = at*(1 - et**2)

#Special case hohmann transfer
nu1 = 0
nu2 = np.pi

nbarm = u.mean_motion(rm, mu_sun)
nbare = u.mean_motion(re, mu_sun)

gam1 = (nu2 - nu1) - nbarm*TOF

print(nbarm)

print("#6 Hohmann Transfer TOF: ", f'{TOF:.4f}' )
print("#7 Initial Phase: ", f'{gam1:.4f}' )

gam2 = gam1 + (nbarm - nbare)*TOF

print("#8 Relative Phase at hohmann transfer out: ", f'{gam2:.4f}' )
print("#9 Hohmann Transfer (back) TOF: ", f'{TOF:.4f}' )

gam3 = -gam2
print("#10 Initial Phase for return trip: ", f'{gam3:.4f}' )

deltaTwait = (gam2 - gam3 + 2*np.pi)/abs(nbarm - nbare)
print("#11 Wait Time on mars before return: ", f'{deltaTwait:.4f}' )

mu_sun = 1.327 * 10**11
re = 1.496 * 10**8 #orbit radius of earth
R_sun = 695700 #actual radius of sun

ve = np.sqrt(mu_sun/re)

print("#12 Earth Circular Velocity: ", f'{ve:.4f}' )

at = (re + R_sun)/2
et = (re - R_sun)/(re + R_sun)
h = np.sqrt(mu_sun*at*(1-et**2))
va = h/re

print("#13 Velocity at Aphelion: ", f'{va:.4f}' )
vinf = ve - va
print("#14 Excess Velocity (relative to earth): ", f'{vinf:.4f}' )

eps_dep = (1/2)*(vinf**2)
a_park = 6378
v0_e = np.sqrt(2*(u.MU_EARTH/a_park + eps_dep))
deltaV = v0_e - 0.465
print("#15 DeltaV", f'{deltaV:.4f}' )




