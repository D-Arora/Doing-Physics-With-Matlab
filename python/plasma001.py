# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 11:07:18 2026

@author: Owner
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def equations_of_motion(t, y, q, m, E, B):
    """
    y: state vector [rx, ry, rz, vx, vy, vz]
    Returns derivatives [vx, vy, vz, ax, ay, az]
    """
    v = y[3:6]
    # Lorentz force: F = q(E + v x B) -> a = F/m
    force = q * (E + np.cross(v, B))
    accel = force / m
    return [v[0], v[1], v[2], accel[0], accel[1], accel[2]]

# --- Parameters ---
q = 1.602e-19       # Charge (Coulombs, e.g., proton)
m = 1.673e-27       # Mass (kg, e.g., proton)
E = np.array([0.0, 1e3, 0.0])   # Electric field (V/m) in y-direction
B = np.array([0.0, 0.0, 0.1])   # Magnetic field (Tesla) in z-direction
E = 0
# Initial conditions: [x, y, z, vx, vy, vz]
r0 = [0.0, 0.0, 0.0]
v0 = [0.0, 1e5, 0.0]
y0 = r0 + v0

# Time span for simulation
t_span = (0, 1e-6)
t_eval = np.linspace(t_span[0], t_span[1], 1200)

# --- Solve ODE ---
sol = solve_ivp(
    equations_of_motion, 
    t_span, 
    y0, 
    args=(q, m, E, B), 
    t_eval=t_eval, 
    rtol=1e-8, 
    atol=1e-10
)

# --- Plotting the Trajectory ---
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection='3d')
ax.plot(sol.y[0], sol.y[1], sol.y[2], label='Trajectory', color='b')
ax.set_title('Charged Particle Motion in E and B Fields')
ax.set_xlabel('X Position (m)')
ax.set_ylabel('Y Position (m)')
ax.set_zlabel('Z Position (m)')
plt.show()