# AOE 3154 
# Homework 5 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u

#2
deltaN = 390
TP = u.timeperiod(deltaN) #days
TP = TP * 24 * u.SIDEREAL_OVER_SOLAR
print("#2 Blue Time Period (Sid. hrs) = ", f'{TP:.4f}')

#5
deltaN = 270
TP = u.timeperiod(deltaN) #days
TP = TP * 24 * u.SIDEREAL_OVER_SOLAR
print("#2 Blue Time Period (Sid. hrs) = ", f'{TP:.4f}')

#9
H = 1032 #km
alpha = 15 #deg
eps = np.arccos((u.R_EARTH+H)*np.sin(15 * np.pi/180)/u.R_EARTH) * 180/np.pi
print("#9 Elevation Angle (deg) = ", f'{eps:.4f}')
beta = 90 - alpha - eps
print("#10 Earth Central Angle (deg) = ", f'{beta:.4f}')
areaPercent = (1-np.cos(beta*np.pi/180))/2 *100
print("#11 Percent Earth visible (%) = ", f'{areaPercent:.4f}')
Tvis = u.timeVisible(u.R_EARTH, H, u.MU_EARTH, beta)/60
print("#12 Time visible (m) = ", f'{Tvis:.4f}')
numSats = 360/(2*beta)
print("#13 Number required = ", f'{numSats:.4f}')

#14
p = 1032 + u.R_EARTH
raan_dot = 360/365.2422 /( 24 * 60 * 60 )#deg/s
mean_motion = np.sqrt(u.MU_EARTH/(p**3)) * 180/np.pi #degrees/s
RHS = -(3/2)*u.J2*((u.R_EARTH/p)**2)*mean_motion
i = np.arccos(raan_dot/RHS) * 180/np.pi
print("#13 Inclination (deg) = ", f'{i:.4f}')
LTAN = (12 - (180-90)/15)
print("#14 LTAN (hours) = ", f'{LTAN:.4f}')





