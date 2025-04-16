# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 10:48:06 2025

@author: alexa
"""


#Import Python files
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from scipy.spatial import cKDTree



pi = np.pi
r = 6371            # radius
x0, y0 = 0, 0     # center of the circle

# Parametric angle
theta = np.linspace(0, 2 * pi, 100)


phi_deg = 30
i = 50  # Inclination of the orbit
theta_rot_deg = 30


phi = np.radians(phi_deg)
i_rad = np.radians(i)
theta_rot = np.radians(theta_rot_deg)

# Define a great circle in 3D with inclination incl
sigma = np.linspace(0, 2 * pi, 500)
x = r * np.cos(sigma)
y = r * np.sin(sigma) * np.cos(i_rad)  # projection squashed by cos(i)

x_orbit = x * np.cos(theta_rot) - y * np.sin(theta_rot)
y_orbit = x * np.sin(theta_rot) + y * np.cos(theta_rot)


if i_rad < pi/2:
    print("prograde")
elif pi/2 < i_rad < pi:
    print("retrograde")
else:
    print("pass")
    
if phi > i_rad:
    print("No launch window")
elif phi == i_rad:
    print("One launch window")
elif phi < i_rad:
    print("Two Lanch windows")
else:
    pass
    
    
phi_comp = pi/2 - phi


# Parametric equations of the circle
x_earth = x0 + r * np.cos(theta)
y_earth = y0 + r * np.sin(theta)

x_var = r * np.sin(phi_comp) * np.cos(theta)
y_var = r * np.sin(phi_comp) * np.sin(theta)

orbit_points = np.vstack((x_orbit, y_orbit)).T
cap_points = np.vstack((x_var, y_var)).T

tree = cKDTree(cap_points)
dists, idxs = tree.query(orbit_points, distance_upper_bound=60)  # tolerance

# Only keep points with valid intersection
intersect_idx = np.where(np.isfinite(dists))[0]
intersects = orbit_points[intersect_idx]



x_array = np.arange(-r,r, 0.1)
x_1 = intersects[4][0]
y_1 = intersects[4][1]
m_grad = y_1/x_1

y_array = x_array*m_grad
norm = np.sqrt(x_array[0]**2 + y_array[0]**2)

x_array = x_array * r / norm
y_array = y_array * r / norm


mask = x_array >= 0
x_array_half = x_array[mask]
y_array_half = y_array[mask]

midpoint = len(x_array) // 2
x_array_half = x_array[:midpoint]
y_array_half = y_array[:midpoint]



print(x_array[0], y_array[0])
# Plot
plt.figure(figsize=(6,6))
plt.plot(x_earth, y_earth, label='Earth')
plt.plot(intersects[:,0], intersects[:,1], 'ko', label='Meridian Crossings')
plt.plot(x_var, y_var, label="Luanch")
plt.plot(x_orbit, y_orbit, 'r--', label=f'Inclined Orbit ({i}°)')
plt.plot(x0, y0, 'ro', label='Center')  # mark the center
plt.plot(x_array_half,y_array_half, 'b')
plt.gca().set_aspect('equal', adjustable='box')  # make it a true circle
plt.grid(True)
#plt.legend()
plt.title("Launch Windows")
plt.xlabel('x')
plt.ylabel('y')
plt.show()











# import cartopy.crs as ccrs
# import matplotlib.pyplot as plt

# # Same as before
# i_deg = 45
# i_rad = np.radians(i_deg)
# lon = np.linspace(-180, 180, 500)
# lat = np.degrees(np.arcsin(np.sin(i_rad) * np.sin(np.radians(lon))))

# # Plot using Cartopy
# fig = plt.figure(figsize=(10, 5))
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.coastlines()
# ax.gridlines(draw_labels=True)

# # Plot the inclined line
# ax.plot(lon, lat, color='red', label=f'{i_deg}° inclined line')
# ax.legend()
# plt.title('Inclined Line on Earth (Plate Carrée Projection)')
# plt.close()



       

