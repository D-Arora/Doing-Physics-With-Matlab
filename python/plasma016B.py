# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:49:14 2026

@author: Owner
"""

import numpy as np
from scipy.integrate import solve_ivp

# --- 1. Constants & Configuration ---
Q = 1.602e-19       # Proton charge (C)
M = 1.673e-27       # Proton mass (kg)
R_COIL = 0.5        # Coil radius (m)
Z_COIL = 1.0        # Coil positions at z = -1.0 and z = +1.0
I_COIL = 5.0e5      # Coil current (A)
MU_0 = 4 * np.pi * 1e-7

def get_B_field(x, y, z):
    """Computes axisymmetric magnetic mirror field using Taylor expansion."""
    def B_axial(z_pos):
        denom = (R_COIL**2 + (z - z_pos)**2)**1.5
        Bz = (MU_0 * I_COIL * R_COIL**2) / (2 * denom)
        dBz_dz = -1.5 * (MU_0 * I_COIL * R_COIL**2 * (z - z_pos)) / (R_COIL**2 + (z - z_pos)**2)**2.5
        return Bz, dBz_dz

    Bz1, dBz_dz1 = B_axial(-Z_COIL)
    Bz2, dBz_dz2 = B_axial(Z_COIL)
    
    Bz = Bz1 + Bz2
    dBz_dz = dBz_dz1 + dBz_dz2
    
    Bx = -0.5 * x * dBz_dz
    By = -0.5 * y * dBz_dz
    return np.array([Bx, By, Bz])

# --- 2. Equations of Motion ---
def lorentz_deriv(t, state):
    x, y, z, vx, vy, vz = state
    v_vec = np.array([vx, vy, vz])
    B_vec = get_B_field(x, y, z)
    a_vec = (Q / M) * np.cross(v_vec, B_vec)
    return [vx, vy, vz, a_vec[0], a_vec[1], a_vec[2]]

# --- 3. Run Simulation ---
init_pos = [0.05, 0.0, 0.0]          # Displaced slightly from axis
init_vel = [1.0e5, 3.0e5, 1.5e5]     # vx, vy, vz (m/s)
initial_state = init_pos + init_vel

t_span = (0, 8e-5)
t_eval = np.linspace(t_span[0], t_span[1], 4000)

sol = solve_ivp(lorentz_deriv, t_span, initial_state, t_eval=t_eval, rtol=1e-7, atol=1e-9)

# --- 4. Dynamic Loss Cone & Pitch Angle Analysis ---
x_vals, y_vals, z_vals = sol.y[0], sol.y[1], sol.y[2]
vx_vals, vy_vals, vz_vals = sol.y[3], sol.y[4], sol.y[5]

# Calculate fields at the center (minimum B) and throat (maximum B) to establish the geometry boundaries
B_mid = np.linalg.norm(get_B_field(0, 0, 0))
B_max = np.linalg.norm(get_B_field(0, 0, Z_COIL))
mirror_ratio = B_max / B_mid

# Theoretical loss cone angle from the midplane (in degrees)
# Particles with a pitch angle smaller than this at the midplane will escape
theta_c_midplane = np.degrees(np.arcsin(np.sqrt(1 / mirror_ratio)))

print(f"--- Magnetic Mirror Properties ---")
print(f"B at center (midplane): {B_mid:.4f} T")
print(f"B at throat (peak):     {B_max:.4f} T")
print(f"Mirror Ratio (Rm):      {mirror_ratio:.2f}")
print(f"Critical Loss Cone Angle (at midplane): {theta_c_midplane:.2f}°\n")

print(f"--- Step-by-Step Trajectory Tracking ---")
# Analyze pitch angle dynamics along the particle's path
for i in range(0, len(t_eval), 500):  # Sample a few points for brevity
    x, y, z = x_vals[i], y_vals[i], z_vals[i]
    v_vec = np.array([vx_vals[i], vy_vals[i], vz_vals[i]])
    v_mag = np.linalg.norm(v_vec)
    
    B_vec = get_B_field(x, y, z)
    B_mag = np.linalg.norm(B_vec)
    
    # Decompose velocity into parallel and perpendicular components relative to local B vector
    v_parallel = np.dot(v_vec, B_vec) / B_mag
    v_perp = np.linalg.norm(v_vec - (v_parallel * B_vec) / B_mag)
    
    # Calculate local pitch angle: angle between velocity vector and magnetic field vector
    pitch_angle = np.degrees(np.arccos(np.abs(v_parallel) / v_mag))
    
    # Calculate local dynamic loss cone angle
    # As B increases toward the throat, the local loss cone widens until it reaches 90 degrees at the mirror point
    if B_mag < B_max:
        local_theta_c = np.degrees(np.arcsin(np.sqrt(B_mag / B_max)))
    else:
        local_theta_c = 90.0

    status = "CONFINED (Reflecting)" if pitch_angle > local_theta_c else "LOST (In Loss Cone)"
    
    print(f"t={t_eval[i]*1e6:5.2f}µs | z={z:5.2f}m | Local B={B_mag:.3f}T | Pitch∠={pitch_angle:5.1f}° | Local Loss Cone ∠={local_theta_c:5.1f}° -> {status}")