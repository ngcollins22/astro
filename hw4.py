# AOE 3154 
# Homework 4 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u

#9-18 or whateva
r = [5000, 4000, 6000] #km
v = [-7, 5, 3] #km

specific_energy = u.specific_energy(r, v)
hvec = np.cross(r, v)

a = -u.MU_EARTH/(2*specific_energy) #km
TP = u.timePeriod(a, u.MU_EARTH) #seconds
mean_motion = ((2*np.pi)/TP) * (60*60*24) /(2*np.pi) #converted to rev/day
print("#11 Mean Motion (rev/day) = ", f'{mean_motion:.4f}')

evec = u.eccenticity_vector(r, v, u.MU_EARTH)
eccentricity = u.mag(evec)
print("#12 Eccentricity = ", f'{eccentricity:.4f}')

nhat = u.unit_line_of_nodes(hvec)
RAAN = u.RAAN(nhat) * 180/np.pi
print("#13 RAAN (degrees) = ", f'{RAAN:.4f}')

inclination = u.inclination(hvec) * 180/np.pi
print("#14 inclination (degrees) = ", f'{inclination:.4f}')

argument_of_periapsis = u.argument_of_periapsis(evec, nhat) * 180/np.pi
print("#15 argument of periapsis (degrees) = ", f'{argument_of_periapsis:.4f}')

p = u.parameter_of_orbit(u.mag(hvec),u.MU_EARTH)
rhere = u.radius_of_orbit(p, eccentricity, 60 * np.pi/180)
print("#16 radius at true anomaly 60 deg (km) = ", f'{rhere:.4f}')

rp = u.radius_of_orbit(p, eccentricity, 0)
altp = rp - u.R_EARTH
print("#17 altitude at periapsis (km) = ", f'{altp:.4f}')

ra = u.radius_of_orbit(p, eccentricity, np.pi)
va = u.mag(hvec)/ra
print("#18 velocity at apogee (km/s) = ", f'{va:.4f}')










