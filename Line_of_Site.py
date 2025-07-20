# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 21:38:49 2025

@author: alexa
"""
#Import Python files
import numpy as np
from numpy import sin, cos, radians


#Data
mu_sun = 1.3271244*10**11
mu_earth = 3.986*10**5
radius_earth = 6378.14
J_2_earth = 1082.63*10**-6
mu_moon = 4902.8
radius_moon = 1737.4
r_mean = 384400


'''
Functions
'''

pi = np.pi


def norm(vec):
    return np.linalg.norm(vec)


def radians(deg):
    rads = deg * (pi/180)
    return rads


def get_mean_anomaly(del_t, M_old, a):
    nu_t = (mu_sun / (a**3))**0.5
    mean_anomaly_t = M_old + nu_t * (del_t)
    
    return mean_anomaly_t

def time_orbit(a, μ):
    T = 2 * np.pi * np.sqrt(a**3 / μ)
    time = T 
    return time

def rotation_matrix(i, Omega, omega):
    cos_O, sin_O = np.cos(Omega), np.sin(Omega)
    cos_i, sin_i = np.cos(i), np.sin(i)
    cos_w, sin_w = np.cos(omega), np.sin(omega)
    
    R = np.array([
        [cos_O * cos_w - sin_O * sin_w * cos_i, -cos_O * sin_w - sin_O * cos_w * cos_i, sin_O * sin_i],
        [sin_O * cos_w + cos_O * sin_w * cos_i, -sin_O * sin_w + cos_O * cos_w * cos_i, -cos_O * sin_i],
        [sin_w * sin_i, cos_w * sin_i, cos_i]
    ])
    
    return R

'''
Newton-Rasphon Method to return the eccentric anomaly to a certain tolerance
'''

#1.1.1 ----------------------------------------------------------------------

def Kepler(e, M, tol = 1e-12, max_i = 1000):
    
    # Guess the solution is similar to theta_M
    E = M
    
    for i in range(max_i):
        
        # Define function in terms of f(E) = 0
        f_E = E - e * np.sin(E) - M
        
        # Define the derivative
        f_prime = 1 - e * np.cos(E)
        
        del_E = f_E / f_prime
        
        # Calculate the new eccentric anomaly
        E_new = E - del_E
        
        # If the error is within some tolerance return theta
        if np.abs(del_E) < tol:
            theta = 2*np.arctan(np.tan(E_new/2) * ((1+e)/(1-e))**(0.5))
            return theta
        
        # Else restart the for loop with a new value
        E = E_new
        


#1.1.3 ----------------------------------------------------------------------

t_0 = 2460705.5*(3600*24)





t_0_days = 2460705.5
days_convert = 3600*24
#1.1.4 ----------------------------------------------------------------------

# Function to return the position and velocity arrays at some time t
def Ephemeris(t, OBJdata, mu):

    time, a, e, i, Omega, omega, mean_anomaly = OBJdata[0:7]
    
    # Find the mean motion
    nu_t = (mu / (a**3))**0.5
    
    # Find the mean anomaly at some time t 
    t = t - t_0_days*days_convert
    mean_anomaly_t = mean_anomaly + nu_t * (t)

    h = np.sqrt(mu * a * (1 - e**2))
    theta_var = Kepler(e, mean_anomaly_t)
    r = a*(1-(e**2))/(1 + e*np.cos(theta_var))
    
    arr_r = np.array([r*np.cos(theta_var), r*np.sin(theta_var), 0])
    arr_v = (mu/h)* np.array([-np.sin(theta_var), e + np.cos(theta_var), 0])
    
    # Perform matrix operations on position and velocity arrays
    R_matrix = rotation_matrix(i, Omega, omega)
    r_ijk = R_matrix @ arr_r
    v_ijk = R_matrix @ arr_v
    
    return r_ijk, v_ijk





t_0_jd = 2460705.5
days_to_sec = 86400
satellite_data = [t_0_jd, 
                  7000,  # semi-major axis in km
                  0.001, 
                  radians(98.7), 
                  radians(120), 
                  radians(87), 
                  radians(45)]  # all angles in radians




def gst_from_jd(jd):
    """
    Compute Greenwich Sidereal Time (GST) in radians at a given Julian Date (jd).
    Returns the angle between the ECI and ECEF frames.
    """
    # Julian centuries since J2000.0
    T = (jd - 2451545.0) / 36525.0

    # GST in seconds (from Vallado)
    theta_sec = 67310.54841 + (876600.0 * 3600 + 8640184.812866) * T \
                + 0.093104 * T**2 - 6.2e-6 * T**3

    # Convert to degrees, then radians
    theta_deg = (theta_sec / 240.0) % 360.0  # 1 sidereal second = 1/240 deg
    theta_rad = np.radians(theta_deg)

    return theta_rad



R_earth = 6378.137  # Earth's equatorial radius in km

def ground_station_eci(lat_deg, lon_deg, alt_m, time_jd):
    lat = radians(lat_deg)
    lon = radians(lon_deg)
    alt = alt_m / 1000.0  # convert to km

    # Convert to ECEF
    x = (R_earth + alt) * cos(lat) * cos(lon)
    y = (R_earth + alt) * cos(lat) * sin(lon)
    z = (R_earth + alt) * sin(lat)
    r_ecef = np.array([x, y, z])

    # Now rotate ECEF to ECI
    theta_gst = gst_from_jd(time_jd)  # in radians
    rotation = np.array([[cos(theta_gst), -sin(theta_gst), 0],
                         [sin(theta_gst),  cos(theta_gst), 0],
                         [0,               0,              1]])
    
    r_eci = rotation @ r_ecef
    return r_eci

def az_el_from_los(LOS, r_gs):
    # Build local ENU basis
    lat = np.arcsin(r_gs[2] / np.linalg.norm(r_gs))
    lon = np.arctan2(r_gs[1], r_gs[0])
    
    sin_lat, cos_lat = np.sin(lat), np.cos(lat)
    sin_lon, cos_lon = np.sin(lon), np.cos(lon)

    # ENU basis vectors
    east  = np.array([-sin_lon, cos_lon, 0])
    north = np.array([-sin_lat*cos_lon, -sin_lat*sin_lon, cos_lat])
    up    = np.array([cos_lat*cos_lon,  cos_lat*sin_lon, sin_lat])

    # Project LOS onto ENU
    e = np.dot(LOS, east)
    n = np.dot(LOS, north)
    u = np.dot(LOS, up)

    # Calculate azimuth and elevation
    azimuth = np.arctan2(e, n) % (2 * np.pi)
    elevation = np.arcsin(u / np.linalg.norm(LOS))
    

    return np.degrees(azimuth), np.degrees(elevation)



t_now = t_0_jd * days_to_sec
r_sat, v_sat = Ephemeris(t_now, satellite_data, mu_earth)

# Example: Auckland NZ
lat_gs = -36.85
lon_gs = 174.76
alt_gs = 40.0  # meters

r_gs = ground_station_eci(lat_gs, lon_gs, alt_gs, t_0_jd)
LOS = r_sat - r_gs
az, el = az_el_from_los(LOS, r_gs)

print(f"Satellite position (ECI): {r_sat}")
print(f"Ground station position (ECI): {r_gs}")
print(f"Line of Sight vector: {LOS}")
print(f"Azimuth: {az:.2f}°, Elevation: {el:.2f}°")