# AOE 3154 
# Homework 8 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u

#1-3
e = 0.6789
M = 0.2024
E0 = M
E1 = u.singleNRKeppler(E0, e, M)
print("#1 Value after first iteration: ", f'{E1:.4f}')

[E, i] = u.solveNRKeppler(M, e, 10**(-4))
print("#2 Find value within tolerance: ", f'{E:.4f}')
print("#3 Iteration Count: ", i)

#4-6
a = 5 #DU
e = 0.4
mu = 1
deltat = 50
mean_motion = u.mean_motion(a, mu)
M = mean_motion*deltat % (2*np.pi)
print("#4 Mean Anomaly (rad): ", f'{M:.4f}')
[E, i] = u.solveNRKeppler(M, e, 1*10**(-4))
print("#5 Eccentric Anomaly (rad): ", f'{E % (2*np.pi):.4f}')
nu = 2*np.pi - u.true_anomaly_Ee(E, e) #this is broken and needs a check in the function
print("#6 True Anomaly (deg): ", f'{nu*180/np.pi:.4f}')


r0vec = [12000,0,0]
v0vec = [1,3,5]
r0 = u.mag(r0vec)
v0 = u.mag(v0vec)

eps = (v0**2)/2 - u.MU_EARTH/r0
a = -u.MU_EARTH/(2*eps)
print("#7 Semimajor Axis: ", f'{a:.4f}')

hvec = np.cross(r0vec, v0vec)
h = u.mag(hvec)
evec = u.eccenticity_vector(r0vec, v0vec, u.MU_EARTH)
e = u.mag(evec)
print("#8 Eccentricity: ", f'{e:.4f}')

nu = u.true_anomaly(evec=evec, rvec=r0vec, vvec=v0vec);
E0 = u.eccentric_anomaly(e, nu);
M0 = u.mean_anomaly(E0, e);

print("#9 E0: ", f'{E0:.4f}')
print("#10 M0: ", f'{M0:.4f}')

deltaT = 6 * 3600 #seconds
n = u.mean_motion(a, u.MU_EARTH);
M = n*deltaT + M0
print("#11 M at 6hrs: ", f'{M % (2*np.pi):.4f}')
[E, i] = u.solveNRKeppler(M, e, 10**(-4))
print("#12 E at 6hrs: ", f'{E % (2*np.pi):.4f}')

dE = E - E0
f = 1 - (a/r0)*(1-np.cos(dE))
g = deltaT - np.sqrt((a**3)/u.MU_EARTH)*(dE - np.sin(dE))

print(f, g)

rvec = f*np.array(r0vec) + g*np.array(v0vec)
r = u.mag(rvec);

fdot = - (np.sqrt(u.MU_EARTH*a)*np.sin(dE))/(r*r0)
gdot = 1 - (a/r)*(1-np.cos(dE))

print("#13-15 R vector at 6 hrs: ", rvec)
print("#16 fdot*1000000", f'{fdot*1000000:.4f}')
print("#17 gdot", f'{gdot:.4f}')

r1 = [3, -2, -1]
r2 = [2, -3, -4]
cross = np.cross(r1, r2)
print(cross)

z = -(3*cross[0] + 2*cross[1])/cross[2]
print("#18 Z = ", f'{z:.4f}')
r3 = [3,2,z]

[p, e, Phat, Qhat, What] = u.gibbs_method(r1,r2,r3);
print("#19 p = ", f'{p:.4f}')
print("#20 e = ", f'{e:.4f}')
print(What)
i = np.arccos(np.dot(What, u.K))
print("#20 i = ", f'{i * 180/np.pi:.4f}')