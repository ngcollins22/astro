# AOE 3154 
# Homework 5 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u

#1-5
mu = 1 #DU3/TU2
r = 2 #DU
i0 = 0 #rad
i1 = np.pi/2 #rad
v = np.sqrt(mu/r)
deltaV = u.planeChange(v, i1-i0) #DU/TU
print("#1 Delta V for plane change: ", f'{deltaV:.4f}')

r1 = 2
r2 = 30
deltaV1 = u.raiseApogee(r1, r2, mu)
print("#2 Radius of apogee: ", f'{r2:.4f}')
print("#3 Delta V to raise apogee: ", f'{deltaV1:.4f}')
v2 = v + deltaV1
ra = r1
rp = r2
va = v2
h = ra*va
vp = h/rp
print("velocity at periapsis: ", f'{vp:.4f}')


deltaV2 = u.planeChange(vp, i1-i0)
print("#4 Delta V to plane change: ", f'{deltaV2:.4f}')

totalDeltaV = 2*np.abs(deltaV1) + np.abs(deltaV2)
print("#5 Total Delta V: ", f'{totalDeltaV:.4f}')

#6-9
rp = 1 #DU periapsis of probe orbit
r = 2 #DU starting circular orbit
# TP are the same -> semimajor axis must be the same
a = r
ra = 2*a - rp #radius of apoapsis of probe orbit
eps1 = u.specific_energy(1,a)
v = u.velocityFromEnergy(eps1, mu, r) #velocity of probe at release
print("#6 Probe velocity: ", f'{v:.4f}')

vp = u.velocityFromEnergy(eps1, mu, rp)
phi = -u.flightPathAngle(vp*rp,r,v) * 180/np.pi
print("#7 Flight path angle at release: ", f'{phi:.4f}')

e2 = 0
p2 = a*(1-e2**2)
et = (ra - rp)/(ra+rp)
pt = a*(1-et**2)

nu = 360 - u.true_anomaly(pt,et,p2,e2)*180/np.pi
#I want the one where flight path angle is negative
print("#8 True Anomaly: ", f'{nu:.4f}')

#the velocities work out to be the same so the maneuver is a plane change through phi
deltaV = np.abs(u.planeChange(v, phi*np.pi/180))
print("#9 DeltaV: ", f'{deltaV:.4f}')

#10-15
e = 0.5
a = 10
E = np.pi/2
M = u.mean_anomaly(E, e)
dA = u.deltaA(a, e, M)
print("#11 Swept Area: ", f'{dA:.4f}')

mu = 10
deltaT = u.TOF(M, mu, a)
print("#12 Time of Flight: ", f'{deltaT:.4f}')

print("#15 Mean Anomaly: ", f'{M:.4f}')