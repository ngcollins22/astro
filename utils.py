import datetime
from scipy import constants
import numpy as np

OMEGA_ECEF = 15 #degrees per hour relative to sun
OMEGA_RAAN = 1.0027379093*OMEGA_ECEF #true degrees per hour
SIDEREAL_OVER_SOLAR = 1.0027379093 #sidereal time/solar time
MU_EARTH = 3.986*(10**5) #gravitational paramter of earth km^3/s^2
R_EARTH = 6378 #km

I = [1,0,0]
J = [0,1,0]
K = [0,0,1]

def arcminutesToDegrees(minutes):
    return minutes/60

def arcsecondsToDegrees(seconds):
    return seconds/3600

def hmsToh(hours, minutes, seconds):
    return hours + arcminutesToDegrees(minutes) + arcsecondsToDegrees(seconds)

def angleTohours(angle):
    return angle/OMEGA_ECEF
    

def dateToHours(month, day, hour, minute):
    date = datetime.datetime(2024, month, day, hour, minute)
    start_of_year = datetime.datetime(2024, 1,1,0,0)
    time_diff = date - start_of_year
    hours_diff = time_diff.total_seconds() / constants.hour
    return hours_diff

def amongus():
    return print(chr(sum(range(ord(min(str(not())))))))

def mag(vec):
    return np.linalg.norm(vec)

#Calculate specific energy from r and v
def specific_energy(r, v):
    return (mag(v)**2)/2 - MU_EARTH/mag(r)

def timePeriod(a, mu):
    return 2*np.pi*np.sqrt((a**3)/mu)

def eccenticity_vector(rvec, vvec, mu):
    hvec = np.cross(rvec, vvec)
    evec = np.cross(vvec, hvec)/mu - rvec/mag(rvec)
    return evec

def unit_line_of_nodes(hvecIJK):
    prod = np.cross(K, hvecIJK)
    return prod/mag(prod)

def RAAN(nhatIJK):
    raanPrime = np.arccos(np.dot(nhatIJK,I))
    if np.dot(nhatIJK,J) >= 0:
        return raanPrime
    else:
        return 2*np.pi - raanPrime
    
def inclination(hvecIJK):
    prod = np.dot(K,hvecIJK)
    return np.arccos(prod/mag(hvecIJK))

def argument_of_periapsis(evecIJK, nhatIJK):
    omega_prime = np.arccos(np.dot(evecIJK,nhatIJK)/mag(evecIJK))
    if np.dot(evecIJK, K) >= 0:
        return omega_prime
    else:
        return 2*np.pi - omega_prime
    
def parameter_of_orbit(h, mu):
    return (h**2)/mu

def radius_of_orbit(p, e, nu):
    return p / (1 + e*np.cos(nu))