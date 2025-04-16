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
theta = np.linspace(0, 2 * pi, 10000)


phi_deg = 28.5 # Launch latitude
i = 50  # Inclination of the orbit
theta_rot_deg = 45 # RAAN


phi = np.radians(phi_deg)
i_rad = np.radians(i)
theta_rot = np.radians(theta_rot_deg)

# Define a great circle in 3D with inclination incl
sigma = np.linspace(0, 2 * pi, 10000)
x = r * np.cos(sigma)
y = r * np.sin(sigma) * np.cos(i_rad)  # projection squashed by cos(i)

x_orbit = x * np.cos(theta_rot) - y * np.sin(theta_rot)
y_orbit = x * np.sin(theta_rot) + y * np.cos(theta_rot)

midpoint2 = len(x_orbit) // 2
x_RAAN = x_orbit[midpoint2]
y_RAAN = y_orbit[midpoint2]

raan_angle = np.arctan2(y_RAAN, x_RAAN) % (2 * np.pi)

x_orbit = x_orbit[midpoint2:]
y_orbit = y_orbit[midpoint2:]


if i_rad < pi/2:
    print("Prograde")
elif pi/2 < i_rad < pi:
    print("Retrograde")
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
dists, idxs = tree.query(orbit_points, distance_upper_bound=1.5)  # tolerance

# Only keep points with valid intersection
intersect_idx = np.where(np.isfinite(dists))[0]
intersects = orbit_points[intersect_idx]


intersect_angles = np.arctan2(intersects[:,1], intersects[:,0]) % (2 * np.pi)
mask_in_range = (intersect_angles <= theta_rot) | (intersect_angles >= pi + theta_rot)
intersects_after = intersects[mask_in_range]


tree_var = cKDTree(intersects_after)
pairs = tree_var.query_pairs(r=5)  # r is the distance threshold

# Find indices to remove (keep only one from each close pair)
to_remove = set()
for i, j in pairs:
    to_remove.add(j)  # or i — just one of them

filtered_points = np.delete(intersects_after, list(to_remove), axis=0)








def split_lines(x_1, y_1, direct):
    x_array = np.arange(-r,r, 0.1)
    
    m_grad = y_1 / x_1
    y_array = x_array*m_grad
    
    norm = np.sqrt(x_array[0]**2 + y_array[0]**2)
    
    midpoint = len(x_array) // 2
    
    x_array = x_array * r / norm
    y_array = y_array * r / norm
    
    if direct == 'up':
        x_array_half = x_array[midpoint:]
        y_array_half = y_array[midpoint:]
    elif direct == 'down':
        x_array_half = x_array[:midpoint]
        y_array_half = y_array[:midpoint]
    
    return x_array_half, y_array_half
        
        
        
    
    
    
# x_1 = filtered_points[0][0]
# y_1 = filtered_points[0][1]
# m_grad = y_1/x_1

# y_array = x_array*m_grad
# norm = np.sqrt(x_array[0]**2 + y_array[0]**2)

# x_array = x_array * r / norm
# y_array = y_array * r / norm


# mask = x_array >= 0
# x_array_half = x_array[mask]
# y_array_half = y_array[mask]


# x_array_half = x_array[:midpoint]
# y_array_half = y_array[:midpoint]


x_array = np.arange(-r,r, 0.1)
x_2 = filtered_points[1][0]
y_2 = filtered_points[1][1] 

m_grad = y_2/x_2

y_array = x_array*m_grad


midpoint = len(x_array) // 2



x_array_down = split_lines(filtered_points[0][0], filtered_points[0][1], 'down')[0]
y_array_down = split_lines(filtered_points[0][0], filtered_points[0][1], 'down')[1]

x_array_up = split_lines(filtered_points[1][0], filtered_points[1][1], 'up')[0]
y_array_up = split_lines(filtered_points[1][0], filtered_points[1][1], 'up')[1]



#print(x_array[0], y_array[0])
# Plot
plt.figure(figsize=(6,6))
plt.plot(x_earth, y_earth, label='Earth')

plt.plot(intersects_after[:,0], intersects_after[:,1], 'ko', label='Meridian Crossings')
plt.plot(x_var, y_var, label="Luanch")
plt.plot(x_orbit, y_orbit, 'r--', label=f'Inclined Orbit ({i}°)')
plt.plot(x0, y0, 'ro', label='Center')  # mark the center
plt.plot(x_array_down,y_array_down, 'b')
plt.plot(x_array_up, y_array_up, 'b')
plt.gca().set_aspect('equal', adjustable='box')  # make it a true circle
plt.grid(True)
#plt.legend()
plt.title("Launch Windows")
plt.xlabel('x')
plt.ylabel('y')
plt.show()






# Calucluations --------------

gamma = np.arcsin(np.cos(i_rad)/np.cos(phi))
print(np.degrees(gamma))

delta = np.arcsin(np.tan(phi)/np.tan(i_rad))

print(np.degrees(delta))



       

