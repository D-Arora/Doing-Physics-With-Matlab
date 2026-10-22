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
import time

tStart = time.time()
plt.close('all')

#%% ALL QUANTITIES GIVEN S.I. UNITS
# Simple magnetic bottle model using two parallel current loops (coils) at z = -d and z = +d
# INPUTS
q = 1.602e-19       # Particle charge (Proton, C)
m = 1.673e-27       # Particle mass (Proton, kg)
Rcoil = 0.5         # Radius of the mirror coils (m)
Zcoil = 1.0         # Distance of coils from origin (m), coils at z = -1.0 and z = +1.0
Icoil = 5.0e5       # Coil current (Amperes)
mu0 = 4*pi*1e-7

tS = 50e-6           # Simulation time (s)       

#%%
# Computes a first-order approximation of an axisymmetric magnetic mirror field.
# Uses the analytical on-axis field Bz(z) and applies a Taylor expansion for Br.
def get_B_field(x, y, z):
    def B_axial(z_pos):
        
        Bm = 1; B0 = 0.4   # Field on axis from a single loop at z_pos
        zm = 0.5
        Bz = B0 #+ (Bm/zm**2) * z**2
        
        #denom = (Rcoil**2 + (z - z_pos)**2)**(1.5)
        #Bz = (mu0 * Icoil * Rcoil**2) / (2 * denom)
        #Bz = Bz*(1+z_pos**2)
        # Derivative dBz/dz for radial field approximation
        #dBz_dz = -1.5 * (mu0 * Icoil * Rcoil**2 * (z - z_pos)) / (Rcoil**2 + (z - z_pos)**2)**(2.5)
        dBz_dz = -2*(Bm/zm**2) * z
        dB_dz = 0
        return Bz, dBz_dz

    # Sum fields from both coils
    Bz1, dBz_dz1 = B_axial(-Zcoil)
    Bz2, dBz_dz2 = B_axial(Zcoil)
    
    Bz = Bz1 + Bz2
    dBz_dz = dBz_dz1 + dBz_dz2
    
    # Linear approximation for radial field from Gauss's Law (div B = 0)
    # Br = -0.5 * r * (dBz/dz)
    Bx = -0.5 * x * dBz_dz
    By = -0.5 * y * dBz_dz
    
    return np.array([Bx, By, Bz])

#%% Solve Equations of Motion EM
# State vector = [x, y, z, vx, vy, vz]
# Returns [vx, vy, vz, ax, ay, az] for solve_ivp.
def EM(t, state):
    x, y, z, vx, vy, vz = state
    v = np.array([vx, vy, vz])
    # Get local magnetic field
    B = get_B_field(x, y, z)
    
    # Acceleration: a = (q/m) * (v x B)
    a = (q / m) * np.cross(v, B)
    
    return [vx, vy, vz, a[0], a[1], a[2]]

# --- 3. Set Initial Conditions & Run Simulation ---
# Start near the center with enough perpendicular velocity to guide gyro-motion, 
# and moderate parallel velocity so it bounces rather than escaping the loss cone.
init_pos = [0.05, 0.0, 0.0]                  # x, y, z (m)
init_vel = [1e5, 1, 1e5]             # vx, vy, vz (m/s)
initial_state = init_pos + init_vel

t_span = (0, tS)                           # Simulate for 80 microseconds
t_eval = np.linspace(t_span[0], t_span[1], 4000) # 4,000 time steps for a smooth orbit

# Solve ODE
sol = solve_ivp(EM, t_span, initial_state, t_eval=t_eval) #, rtol=1e-7, atol=1e-9)

# Extract coordinates
t = t_eval
x, y, z = sol.y[0], sol.y[1], sol.y[2]
vx, vy, vz = sol.y[3], sol.y[4], sol.y[5]
v = (vx**2+vy**2+vz**2)**0.5
vR = (vx**2+vy**2)**0.5

#%%  GRAPHICS
plt.rcParams['font.size'] = 10
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (4.3,4))
ax.set_xlabel(' x  [ m ]')
ax.set_ylabel('y  [ m ]')
ax.grid()
ax.plot(x[0],y[0],'go',ms = 8)
ax.plot(x[-1],y[-1],'ro',ms = 8)
ax.plot(x,y,'b',lw = 2)
ax.axis("equal")
fig1.tight_layout()

plt.rcParams['font.size'] = 10
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (4.3,4))
ax.set_xlabel(' t  [ m ]')
ax.set_ylabel('x,y  [ m ]')
ax.grid()
ax.plot(t,x,'b',lw = 2)
ax.plot(t,y,'r',lw = 2)
fig2.tight_layout()

plt.rcParams['font.size'] = 10
fig3, ax = plt.subplots(nrows=1, ncols=1,figsize = (4.3,4))
ax.set_xlabel(' t  [ m ]')
ax.set_ylabel('z  [ m ]')
ax.grid()
ax.plot(t,z,'b',lw = 2)
fig3.tight_layout()

plt.rcParams['font.size'] = 10
fig4 = plt.figure(figsize=(5.5,4))
# 3D trajectory
ax = fig4.add_subplot(111, projection="3d")
ax.plot(x, y, z, color="blue")
ax.plot(x[0],y[0],z[0],'go',ms = 6)
ax.plot(x[-1],y[-1],z[-1],'ro',ms = 6)         
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
fig4.tight_layout()

plt.rcParams['font.size'] = 12
fig5, ax = plt.subplots(nrows=1, ncols=1,figsize = (6,4))
ax.set_xlabel(' t  [ s ]')
ax.set_ylabel('v  [ m/s ]')
ax.grid()
ax.plot(t,vz,lw = 2,label = 'vz')
ax.plot(t,v,lw = 2,label = 'v')
ax.plot(t,vR,lw = 2,label = 'vR')


ax.legend()
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
fig5.tight_layout()


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)

#%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png') 
fig3.savefig('a3.png') 
fig4.savefig('a4.png')
fig5.savefig('a5.png')  
'''
