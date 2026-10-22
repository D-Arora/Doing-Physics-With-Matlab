# -*- coding: utf-8 -*-
"""
Oct 2026         plasma016py

[3D] Motion of charged particles in magnetic and electric fields
  
Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/plasma001.py
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, tan, sin 



#%% ALL QUANTITIES GIVEN S.I. UNITS
# Simple magnetic bottle model using two parallel current loops (coils) at z = -d and z = +d
# INPUTS
Q = 1.602e-19       # Particle charge (e.g., Proton, C)
M = 1.673e-27       # Particle mass (e.g., Proton, kg)
R_COIL = 0.5        # Radius of the mirror coils (m)
Z_COIL = 1.0        # Distance of coils from origin (m), coils at z = -1.0 and z = +1.0
I_COIL = 5.0e5      # Coil current (Amperes)
MU_0 = 4 * np.pi * 1e-7

def get_B_field(x, y, z):
    """
    Computes a first-order approximation of an axisymmetric magnetic mirror field.
    Uses the analytical on-axis field Bz(z) and applies a Taylor expansion for Br.
    """
    def B_axial(z_pos):
        # Field on axis from a single loop at z_pos
        denom = (R_COIL**2 + (z - z_pos)**2)**(1.5)
        Bz = (MU_0 * I_COIL * R_COIL**2) / (2 * denom)
        # Derivative dBz/dz for radial field approximation
        dBz_dz = -1.5 * (MU_0 * I_COIL * R_COIL**2 * (z - z_pos)) / (R_COIL**2 + (z - z_pos)**2)**(2.5)
        return Bz, dBz_dz

    # Sum fields from both coils
    Bz1, dBz_dz1 = B_axial(-Z_COIL)
    Bz2, dBz_dz2 = B_axial(Z_COIL)
    
    Bz = Bz1 + Bz2
    dBz_dz = dBz_dz1 + dBz_dz2
    
    # Linear approximation for radial field from Gauss's Law (div B = 0)
    # Br = -0.5 * r * (dBz/dz)
    Bx = -0.5 * x * dBz_dz
    By = -0.5 * y * dBz_dz
    
    return np.array([Bx, By, Bz])

# --- 2. Define Equations of Motion ---
def lorentz_deriv(t, state):
    """
    Returns [vx, vy, vz, ax, ay, az] for solve_ivp.
    state vector = [x, y, z, vx, vy, vz]
    """
    x, y, z, vx, vy, vz = state
    v_vec = np.array([vx, vy, vz])
    
    # Get local magnetic field
    B_vec = get_B_field(x, y, z)
    
    # Acceleration: a = (q/m) * (v x B)
    a_vec = (Q / M) * np.cross(v_vec, B_vec)
    
    return [vx, vy, vz, a_vec[0], a_vec[1], a_vec[2]]

# --- 3. Set Initial Conditions & Run Simulation ---
# Start near the center with enough perpendicular velocity to guide gyro-motion, 
# and moderate parallel velocity so it bounces rather than escaping the loss cone.
init_pos = [0.05, 0.0, 0.0]                  # x, y, z (m)
init_vel = [1.0e5, 3.0e5, 1.5e5]             # vx, vy, vz (m/s)
initial_state = init_pos + init_vel

t_span = (0, 8e-5)                           # Simulate for 80 microseconds
t_eval = np.linspace(t_span[0], t_span[1], 4000) # 4,000 time steps for a smooth orbit

# Solve ODE
sol = solve_ivp(lorentz_deriv, t_span, initial_state, t_eval=t_eval, rtol=1e-7, atol=1e-9)

# Extract coordinates
x_vals, y_vals, z_vals = sol.y[0], sol.y[1], sol.y[2]

plt.plot(t_eval,x_vals)


