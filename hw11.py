# AOE 3154 
# Homework 10 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u


#2

rv = 0.723 #AU
mu_sun = 39.4784 #AU3/yr2

vesc = np.sqrt(2*mu_sun/rv)

print("#2 Escape Velocity: ", f'{vesc:.2f}')

vcirc_venus = np.sqrt(mu_sun/rv)
vinf = vesc - vcirc_venus
print("#3 V infinity: ", f'{vinf:.2f}')

vinf = vinf * 4.744
vcirc_venus = vcirc_venus * 4.744 # into km/s
mu_venus = 3.249 * 10**5
adep = -mu_venus/(vinf**2)
print("#4 Departure Hyperbola SMA: ", f'{adep:.2f}')

Rv = 6052
Hmin = 400
rpmin = Rv + Hmin
ehyp = -rpmin/adep + 1
delta = 2 * np.arcsin(1/ehyp)
print("#5 Turning angle (deg): ", f'{delta*180/np.pi:.2f}')

b = -adep*np.sqrt(ehyp**2 - 1)
print("#6 Impact parameter (km): ", f'{b:.2f}')

deltaVGA = 2*vinf*np.sin(delta/2)/4.744
vout = vesc
gam = (np.pi - delta)/2
beta = np.pi - gam

deltaphi = np.arcsin(deltaVGA/vout * np.sin(beta))
alpha = np.pi - deltaphi - beta
vin_sun = np.sin(alpha)/np.sin(beta) * vout

print("#7 GA Delta V (AU/yr): ", f'{deltaVGA:.2f}')
print("#8 Incoming Heliocentric Velocity (AU/yr): ", f'{vin_sun:.2f}')


