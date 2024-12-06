import datetime
import scipy.constants as constants
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import webbrowser

OMEGA_ECEF = 15 #degrees per hour relative to sun
OMEGA_RAAN = 1.0027379093*OMEGA_ECEF #true degrees per hour
SIDEREAL_OVER_SOLAR = 1.0027379093 #sidereal time/solar time
MU_EARTH = 3.986*(10**5) #gravitational paramter of earth km^3/s^2
R_EARTH = 6378 #km
M_EARTH = 5.972*(10**24) #kg
G = 6.67430*(10**-20) # SI units
M_LUNA = 7.34767309*(10**22) #kg
J2 = 1.082*(10**-3) # dimensionless

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

def hoursToAngle(hours):
    return hours*OMEGA_ECEF
    

def dateToHours(month, day, hour, minute):
    date = datetime.datetime(2024, month, day, hour, minute)
    start_of_year = datetime.datetime(2024, 1,1,0,0)
    time_diff = date - start_of_year
    hours_diff = time_diff.total_seconds() / constants.hour
    return hours_diff

def amongus():
    return print(chr(sum(range(ord(min(str(not())))))))

# Calculate magnitude of vector
def mag(vec):
    return np.linalg.norm(vec)

#Calculate specific energy from r and v
def specific_energy(r, v, mu):
    return (mag(v)**2)/2 - mu/mag(r)

def specific_energy(mu, a):
    return -mu/(2*a)


#Calculate Time period from a and mu
def timePeriod(a, mu):
    return 2*np.pi*np.sqrt((a**3)/mu)

#Caclulate eccentricity vector from r, v and mu
def eccenticity_vector(rvec, vvec, mu):
    hvec = np.cross(rvec, vvec)
    evec = np.cross(vvec, hvec)/mu - rvec/mag(rvec)
    return evec

#Calculate nhat (unit line of nodes) from angular momentum vector
def unit_line_of_nodes(hvecIJK):
    prod = np.cross(K, hvecIJK)
    return prod/mag(prod)

#Calculate RAAN from nhat
def RAAN(nhatIJK):
    raanPrime = np.arccos(np.dot(nhatIJK,I))
    if np.dot(nhatIJK,J) >= 0:
        return raanPrime
    else:
        return 2*np.pi - raanPrime
    
#Calculate inclination from angular momentum vector
def inclination(hvecIJK):
    prod = np.dot(K,hvecIJK)
    return np.arccos(prod/mag(hvecIJK))

#Calculate argument of periapsis from eccentricity vector and nhat
def argument_of_periapsis(evecIJK, nhatIJK):
    omega_prime = np.arccos(np.dot(evecIJK,nhatIJK)/mag(evecIJK))
    if np.dot(evecIJK, K) >= 0:
        return omega_prime
    else:
        return 2*np.pi - omega_prime

#Calculate parameter of orbit from h and mu   
def parameter_of_orbit(h, mu):
    return (h**2)/mu

#Calculate current radius of orbit from p, e, and true anomaly nu
def radius_of_orbit(p, e, nu):
    return p / (1 + e*np.cos(nu))

def true_anomaly(evec, rvec, vvec):
    ehat = normalize(evec)
    rhat = normalize(rvec)
    nuprime = np.arccos(np.dot(ehat, rhat))
    if(np.dot(rvec, vvec) < 0):
        return 2*np.pi - nuprime
    else:
        return nuprime

def normalize(vec):
    return vec/mag(vec)

def mean_motion(TP): #rev/day
    return ((2*np.pi)/TP) * (60*60*24) /(2*np.pi)

#Calculate orbital elements from r, v vectors and mu
def orbital_elements(rvec, vvec, mu):
    print("Disregard outputs if any conditions are violated")
    evec = eccenticity_vector(rvec, vvec, mu)
    e = mag(evec) #need to check if this is <1, >1, 0, etc.
    hvec = np.cross(rvec, vvec)
    nhat = unit_line_of_nodes(hvec)
    i = inclination(hvec)
    if(i == 0 and e == 0):
        print("Orbit is circiular and equitorial")
    elif(e == 0):
        print("Orbit is circular")
    elif(i == 0):
        print("Orbit is equatorial")
    else: #Can calculate orbital elements as normal
        raan = RAAN(nhat)
        omega = argument_of_periapsis(evec,nhat)
        eps = specific_energy(mag(rvec),mag(vvec))
        a = -mu/(2*eps)
        nu = true_anomaly(evec, rvec, vvec)
        return [a, e, i, raan, omega, nu]


def timeperiod(deltaN):
    return (1/SIDEREAL_OVER_SOLAR)*(1-deltaN/360)

def timeVisible(R, H, mu, beta):
    return 2 * np.pi * np.sqrt((R+H)**3 /mu)*(2*beta)/360

def planeChange(v, deltaTheta):
    return 2 * v * np.sin(deltaTheta/2)

def raiseApogee(r1, r2, mu):
    return np.sqrt(2*mu*(1/r1 - 1/(r1+r2))) - np.sqrt(mu/r1)

def velocityFromEnergy(eps, mu, r):
    return np.sqrt(2*(eps + mu/r))

def flightPathAngle(h, r, v):
    return np.arccos(h/(r*v))

def true_anomaly_pe(pt, et, p2, e2):
    return np.arccos((pt-p2)/(et*p2 - e2*pt))

def mean_anomaly(E, e):
    return E - e*np.sin(E)

def deltaA(a, e, M):
    return 0.5*np.sqrt(1 - e**2)*(a**2)*M

def TOF(M, mu, a):
    return M/np.sqrt(mu/(a**3))

def singleNRKeppler(Ei, e, M):
    Mi = mean_anomaly(Ei, e)
    Ei1 = Ei - (Mi-M)/(1 - e*np.cos(Ei))
    return Ei1

def solveNRKeppler(M, e, tol):
    Ei = M
    Mi = mean_anomaly(Ei, e)
    i = 0
    while(np.abs(Mi - M) > tol):
        Ei = singleNRKeppler(Ei, e, M)
        Mi = mean_anomaly(Ei, e)
        i = i + 1
    return [Ei, i]

def mean_motion(a, mu):
    return np.sqrt(mu/(a**3))

def true_anomaly_Ee(E, e):
    cosNu = (np.cos(E) - e)/(1 - e*np.cos(E))
    return np.arccos(cosNu)


def eccentric_anomaly(e, nu):
    return np.arccos((e + np.cos(nu))/(1 + e*np.cos(nu)))

def lagrange_coeffs_2D(r0, v0, r, v):
    h = mag(np.cross(r0, v0))
    f = (r[0]*v0[1] - r[1]*v0[0])/h
    g = (r[1]*r0[0] - r[0]*r0[1])/h
    fdot = (v[0]*v0[1] - v[1]*v0[0])/h
    gdot = (r0[0]*v[1] - r0[1]*v[1])/h
    return [f, g, fdot, gdot]

def gibbs_method(r1, r2, r3):
    Dvec = np.cross(r2,r3) + np.cross(r3, r1) + np.cross(r1,r2)
    Nvec = mag(r1)*np.cross(r2,r3) + mag(r2)*np.cross(r3,r1)+mag(r3)*np.cross(r1,r2)
    Svec = (mag(r2)-mag(r3))*np.array(r1) + (mag(r3)-mag(r1))*np.array(r2) + (mag(r1)-mag(r2))*np.array(r3)
    D = mag(Dvec)
    N = mag(Nvec)
    S = mag(Svec)
    p = N/D
    e = S/D
    Qhat = Svec/S
    What = Nvec/N
    Phat = np.cross(Qhat,What)
    return [p, e, Phat, Qhat, What]

def launch_azimuth(Ls, i):
    if Ls >= 0:
        return np.arcsin(np.cos(i)/np.cos(Ls))
    else:
        return np.pi - np.arcsin(np.cos(i)/np.cos(Ls))

def launch_location_angle(beta, i, Ls):
    if Ls >= 0:
        return np.arccos(np.cos(beta)/np.sin(i))
    else:
        return -np.arccos(np.cos(beta)/np.sin(i))
    
def apparent_velocity(vin, veast, beta):
    product = vin**2 + veast**2 - 2*vin*veast*np.cos(np.pi/2 - beta)
    return np.sqrt(product)

def apparent_launch_azimuth(beta, vin, vapp):
    return np.arcsin((vin/vapp)*np.sin(np.pi/2 - beta)) - np.pi/2

def hohmann_rendevous(achaser, atarget, mu, phi_init):
    print("Hohmann Rendevous from a =", achaser, "to a =", atarget, "** WATCH UNITS **")
    at = (achaser + atarget)/2
    print(" at =", at)
    TOF = np.pi*np.sqrt(at**3 / mu)
    print(" TOF b/w burns =", f'{TOF:.4f}')
    mean_motion_target = np.sqrt(mu/(atarget**3))
    alpha_lead = mean_motion_target*TOF
    print(" alpha lead (rad) =", f'{alpha_lead:.4f}')
    phi_imp = np.pi - alpha_lead
    print(" phase at impulse (rad) =", f'{phi_imp:.4f}')
    mean_motion_chaser =  np.sqrt(mu/(achaser**3))
    k = 0
    deltaT_wait = (phi_imp - phi_init + 2*k*np.pi)/(mean_motion_target - mean_motion_chaser)
    print(" wait time = ", f'{deltaT_wait:.4f}')
    return [TOF, alpha_lead, phi_imp, deltaT_wait]

def p_iteration_method(r1v, r2v, TOF, eps, mu, way):
    # Initial constants
    r1 = mag(r1v)
    r2 = mag(r2v)
    Nu = np.arccos(np.dot(r1v, r2v)/(r1*r2))
    
    if way == 1 :
        Nu = 2*np.pi - Nu
    cNu = np.cos(Nu)
    k = r1*r2*(1-cNu)
    m = r1*r2*(1+cNu)
    l = r1 + r2

    # Initial p-guess
    pi = k/(l + np.sqrt(2*m))
    pii = k/(l - np.sqrt(2*m))
    p0 = (pi + pii)/2
    pn = p0
    pn1 = p0
    TOFn = 0
    while(np.abs(TOFn - TOF) >= eps):
        pn = pn1
        deltaEn = 2 * np.arccos((pn*l - k)/(2*pn*np.sqrt(r1*r2)*np.cos(Nu/2)))
        an = (m*k*pn)/((pn**2)*(2*m - l**2) + 2*k*l*pn - k**2)
        gn = (r1*r2*np.sin(Nu))/np.sqrt(mu*pn)
        TOFn = gn + np.sqrt((an**3)/mu)*(deltaEn - np.sin(deltaEn))
        #print(TOFn)
        #print(pn)

        dTOFdp = -gn/(2*pn) - (3/2)*an*(TOFn - gn)*((k**2 + (2*m - l**2)*(pn**2))/(m*k*(pn**2))) + np.sqrt((an**3)/mu)*(2*k*np.sin(deltaEn))/(pn*(k - l*pn))
        pn1 = pn + (TOF - TOFn)/dTOFdp
    # done!!!!
    #print(TOFn)

    f = 1 - (r2/pn)*(1-cNu)
    g = (r1*r2*np.sin(Nu))/np.sqrt(mu*pn)
    fdot = np.sqrt(mu/pn)*np.tan(Nu/2)*((1-cNu)/pn - 1/r2 - 1/r1)
    gdot = 1 - (r1/pn)*(1-cNu)
    v1v = (np.array(r2v) - f*np.array(r1v))/g
    v2v = fdot*np.array(r1v) + np.array(v1v)*gdot

    return [pn, v1v, v2v]


def RAAN_dot(p, mean_motion, i):
    return (-3/2)*J2*((R_EARTH/p)**2)*mean_motion*np.cos(i)


def walker(t,p,f, star):
    if(star == 1):
        deltaRAAN = (np.pi)/p
    else:
        deltaRAAN = (2*np.pi)/p
    
    satsperplane = t/p
    deltaNuInOrbit = (2*np.pi)/satsperplane
    if(star == 1):
        deltaNuPhase = (np.pi)/t * f
    else:
        deltaNuPhase = (2*np.pi)/t * f

    return [deltaRAAN, deltaNuInOrbit, deltaNuPhase]

def rsoi(am, m, M):
    return am*((m/M)**(2/5))

def chi(m1, m2, m3):
    p = [m1+m2,3*m1+2*m2,3*m1+m2,-(m2+3*m3),-(2*m2+3*m3),-(m2+m3)]
    roots = np.roots(p)
    for root in roots:
        if ~np.iscomplex(root):
            return np.real(root)
        
def find_third_mass(m1, m2, chi):
    den = (3*chi**2 + 3*chi + 1)
    num = (m1+m2)*(chi**5) + (3*m1 + 2*m2)*(chi**4) + (3*m1 + m2)*(chi**3) - m2*(chi**2) - 2*m2*chi - m2
    return num/den
    