# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 14:51:57 2026

@author: Owner
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, tan, sin 

# Define constants
q = 1.0  # charge
m = 1.0  # mass


# Define a simple magnetic mirror field B(z) along the z-axis
def get_B(z, B0=1.0, B_ratio=3.0, L=5.0):
  # Mirror profile: weaker in the middle (z=0), stronger at the ends
  Bz = B0 * (1.0 + (B_ratio - 1.0) * (z / L) ** 2)
  # Simple radial component from div B = 0 approximation for small r
  return np.array([0.0, 0.0, Bz])


# Lorentz force function: dv/dt = (q/m) * (v x B)
def equations_of_motion(t, state):
  # state = [x, y, z, vx, vy, vz]
  x, y, z, vx, vy, vz = state
  r = np.array([x, y, z])
  v = np.array([vx, vy, vz])

  # Magnetic field at position z (neglecting off-axis complexity for basic 1D mirror proxy)
  # For a proper off-axis field, B depends on x, y, z
  B = get_B(z)

  # Lorentz force F = q (v x B)
  F = q * np.cross(v, B)
  a = F / m

  return [vx, vy, vz, a[0], a[1], a[2]]


# Initial conditions: [x, y, z, vx, vy, vz]
# Give it some parallel (vz) and perpendicular (vxy) velocity
v_perp = 0.5
v_parallel = 0.8
state0 = [0.0, 0.0, -4.0, v_perp, 0.0, v_parallel]

t_span = (0, 30)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

sol = solve_ivp(
    equations_of_motion, t_span, state0, t_eval=t_eval, method='RK45'
)

# Plot trajectory (z vs x)
plt.figure(figsize=(8, 4))
plt.plot(sol.y[2], sol.y[0], label='Particle Trajectory')
plt.xlabel('Axial Position z')
plt.ylabel('Radial Position x')
plt.title('Charged Particle Trapping in a Magnetic Mirror')
plt.axvline(-5.0, color='r', linestyle='--', alpha=0.5, label='Mirror Throat')
plt.axvline(5.0, color='r', linestyle='--', alpha=0.5)
plt.grid(True)
plt.legend()
plt.show()



#%%
L = 5
z = linspace(-L,L,200)
B0 = 1; BR = 2
BZ = B0 * (1.0 + 2 * (z / L) ** 2)
plt.plot(z,BZ)