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

vin_sun = np.sqrt(vout**2 + deltaVGA**2 - 2*vout*deltaVGA*np.cos(gam))

print("#7 GA Delta V (AU/yr): ", f'{deltaVGA:.2f}')
print("#8 Incoming Heliocentric Velocity (AU/yr): ", f'{vin_sun:.4f}')

#11-14

m1 = 68
m2 = 47
m3 = 9

chi = u.chi(m1,m2,m3)
print("#11 Chi-Value: ", f'{chi:.2f}')

m1 = 26
m2 = 7
m3 = u.find_third_mass(m1, m2, chi)
print("#12 Third Mass:", f'{m3:.2f}')

r13 = 51
r12 = r13/(1 + chi)

print("#13 R 1->2:", f'{r12:.2f}')

r1 = - (m2 * r12 + m3*r13)/(m1 + m2 + m3)
