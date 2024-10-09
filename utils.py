import datetime
from scipy import constants
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
G = 6.67430*10**(-11) # SI units

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
def specific_energy(r, v):
    return (mag(v)**2)/2 - MU_EARTH/mag(r)


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
    


def visualize(p, e):

    webbrowser.register('chrome', None, webbrowser.GenericBrowser('chrome'))
    pio.renderers.default = 'browser'

    nu = np.linspace(0, 2*np.pi, 360)
    r = radius_of_orbit(p, e, nu)

    fig = go.Figure()

    fig.add_trace(go.Scatterpolar(
        r=r,  # radius values
        theta=nu,  # angle values in degrees
        mode='lines',  # Display as lines
        name='Polar Plot'
    ))

    # Customize layout
    fig.update_layout(
        polar=dict(
            radialaxis=dict(visible=True),
            angularaxis=dict(visible=True)
        ),
        showlegend=True
    )

    # Show the plot
    fig.show()