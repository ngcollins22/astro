# AOE 3154 
# Homework 8 - Nathan Collins

import numpy as np
import scipy as sp
import utils as u

#1
Ls = 34.6320;
mini = Ls;
maxi = 180 - Ls;

print("#1 Min inclination: ", f'{mini:.4f}')
print("#2 Max inclination: ", f'{maxi:.4f}')

veast = 465*np.cos(Ls * np.pi/180)/1000;

print("#3 Eastward Velocity: ", f'{veast:.4f}')

i = 97.7;
H = 567; #km

betaAN = u.launch_azimuth(Ls * np.pi/180, i * np.pi/180)*(180/np.pi)
print("#4 Launch Azimuth AN: ", f'{betaAN:.4f}')

deltaAN = u.launch_location_angle(betaAN*np.pi/180, i * np.pi/180, Ls*np.pi/180)*180/np.pi
print(deltaAN)
RAAN = 105
print("#5 LWST (AN) (hours): ", f'{(RAAN - deltaAN)/15:.4f}')

a = u.R_EARTH + H #km
vcirc = np.sqrt(u.MU_EARTH/a) #km/s
print("#6 Circular/Intertial Velocity: ", f'{vcirc:.4f}')

vapp = u.apparent_velocity(vcirc, veast, betaAN*(np.pi/180))
print("#7 Apparent Launch Velocity: ", f'{vapp:.4f}')

betaANapparent = u.apparent_launch_azimuth(betaAN*(np.pi/180), vcirc, vapp) * 180/np.pi
print("#8 Apparent Launch Azimuth: ", f'{betaANapparent:.4f}')

# 10 - 15

hchaser = 350 #km
htarget = 420 #km
achaser = u.R_EARTH + hchaser
atarget = u.R_EARTH + htarget

phi_init = np.pi
[TOF, _, phi_imp, waitTime] = u.hohmann_rendevous(achaser, atarget, u.MU_EARTH, phi_init)
print("#10 TOF b/w burns (min): ", f'{TOF/60:.4f}')
print("#11 phi impulse (rad): ", f'{phi_imp:.4f}')
print("#12 wait time to initiate (hr): ", f'{waitTime/3600:.4f}')

k = 1
phi_init = 2.9
aphase = atarget*((1 - phi_init/(2*k*np.pi))**(2/3))
print("#13 phasing orbit semimajor axis (km): ", f'{aphase:.4f}')

timeCap = 20 * 3600 #s
mean_motion_target = np.sqrt(u.MU_EARTH/(atarget**3))

klimit = (phi_init + timeCap*mean_motion_target)/(2*np.pi)
print("#14 maximum number of phasing orbits: ", np.floor(klimit))

aphase = atarget*((1 - phi_init/(2* np.floor(klimit)*np.pi))**(2/3))
vt = np.sqrt(u.MU_EARTH*(2/atarget - 1/aphase))
vcirc = np.sqrt(u.MU_EARTH/atarget)
deltaV = (vcirc - vt)*2
print("#15 total deltaV required (km/s): ", f'{deltaV:.4f}')

r1v = [7000, 0, 0]
v1v0 = [0, 5.336, 5.336]
r2v = [-4370, -3235, -4081]
v2v0 = [5.81, -3.242, -3.882]

[p, v1v, v2v] = u.p_iteration_method(r1v, r2v, 3600, 10**(-3), u.MU_EARTH, 1)
print("#16 connecting orbit parameter (km): ", f'{p:.4f}')

deltaV1v = v1v0 - v1v
deltaV2v = v2v0 - v2v
deltaV = u.mag(deltaV1v) + u.mag(deltaV2v)
print("#17 Total deltaV (km/s): ", f'{deltaV:.4f}')

