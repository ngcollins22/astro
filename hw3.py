# AOE 3154 
# Homework 3 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u

# 1
phiH = 19 + u.arcminutesToDegrees(49) + u.arcsecondsToDegrees(20) # in degrees
lHW = 155 + u.arcminutesToDegrees(28) + u.arcsecondsToDegrees(30) # in degrees

dLMT = -lHW/u.OMEGA_ECEF #EAST, in hours
print("#1 Offset from Local Time (hours) ", f'{dLMT:.3f}')

#2
altairRA = u.hmsToh(19, 50, 47) # in hours
thetaG = (altairRA - dLMT) % 24 # in hours
print("#2 Greenwich Sidereal Time (hours) ", f'{thetaG:.3f}')

#3
thetaG0 = u.hmsToh(6,41,33.918) # in hours
hoursSince = u.dateToHours(9, 20, 0, 0) # in hours
thetaG1 = (thetaG0 + u.SIDEREAL_OVER_SOLAR*hoursSince) % 24 # in hours
print("#3 Greenwich Sidereal Time ", f'{thetaG1:.3f}')

#4
thetaLST = thetaG
deltaT = (thetaLST - thetaG1 - u.angleTohours(phiH))/u.SIDEREAL_OVER_SOLAR
hoursDeltaT = deltaT - (deltaT % 1)
minutesDeltaT = (deltaT % 1) * 60
print("#4 UTC Time on Sep 20 hours: ", f'{hoursDeltaT:.3f}')
print("#5 UTC Time on Sep 20 minutes: ", f'{minutesDeltaT:.3f}')

haSThrs = (hoursDeltaT - 10) + 24
print("#6 Hawaiiain-Aleutian Standard Time ", f'{haSThrs:.0f}', ":", f'{minutesDeltaT:.0f}')



#15, 16, 17
mu = 3.986 * 10**5 #km3/s2
p = 6000 #km
e = 0.5
nu = 30 #degrees

v_p = np.sqrt(mu/p)*(-np.sin(nu * np.pi/180))
v_q = np.sqrt(mu/p)*(e + np.cos(nu * np.pi/180))

print("#15 P component: ", f'{v_p:.3f}')
print("#16 Q component: ", f'{v_q:.3f}')


print(chr(sum(range(ord(min(str(not())))))))