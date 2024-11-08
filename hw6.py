# AOE 3154 
# Homework 5 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u

#1
msat = 10 #kg
Rsat = 0.25 #m
Cd=2.2
H = 400 #km
rho = 4 * 10**-12 #kg/m^3
print(rho)
a = u.R_EARTH + H

v=np.sqrt(u.MU_EARTH/a)*1000 #m/s
Aram = np.pi*(Rsat**2)
Fd = (1/2)*rho*Cd*Aram*(v**2) #N
ad = Fd/msat * 1000000 #microm/s^2
print("#1 Acceleration due to Drag (micrometer/s^2)=", f'{ad:.4f}')

eps_dot = -Fd*v/msat * 1000**-2 #km^2/s^3
delta_eps = 24*60*60 * eps_dot
print("#2 Change in specific energy (km^2/s^2)=", f'{delta_eps:.4f}')

eps0 = -u.MU_EARTH/(2*a)
eps1 = eps0 + delta_eps
r1 = -u.MU_EARTH/(2*eps1) 
deltaR = a - r1
print("#3 Change in radius km =", f'{deltaR:.4f}')

mfinal = 50
mfuel = 20
m0 = mfinal + mfuel
Isp = 212
g0 = 9.81
deltaV = Isp * g0 * np.log(m0/mfinal)
print("#4 DeltaV (m/s) =", f'{deltaV:.4f}')

deltaV = 800
massFraction = np.exp(deltaV/(Isp*g0))
print("#5 Mass fraction =", f'{massFraction:.4f}')

mDry2 = 51
mWet2 = massFraction * mDry2

mFuel2 = mWet2 - mDry2

extraFuel = mFuel2-mfuel
print("#6 Additional Fuel =", f'{extraFuel:.4f}')


r = 7000
ra = 7100
rp = 7000

a_before = (ra + rp)/2
eps_before = -u.MU_EARTH/(2*a_before)
vp = np.sqrt(u.MU_EARTH*(2/rp - 1/r))
v_before = np.sqrt(2*(eps_before + u.MU_EARTH/rp))
dv = vp - v_before
print("#8 DeltaV (m/s) =", f'{dv*1000:.4f}')


mu = 1
R1 = 28
R2 = 2

dv1 = np.sqrt(2*mu*(1/R1 - 1/(R1+R2))) - np.sqrt(mu/R1)
dv2 = np.sqrt(mu/R2) - np.sqrt(2*mu*(1/R2 - 1/(R1+R2)))


print("#12 DeltaV 1 (DU/TU) =", f'{dv1:.4f}')
print("#13 DeltaV 2 (DU/TU) =", f'{dv2:.4f}')

at = (R1 + R2)/2


TOF = np.pi * np.sqrt(at**3 / mu)

print("#14 TOF (TU)", f'{TOF:.4f}')






