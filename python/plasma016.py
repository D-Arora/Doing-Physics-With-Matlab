# -*- coding: utf-8 -*-
"""
Oct 2026         plasma016.py

[3D] Motion of charged particles in magnetic and electric fields
Magnetic mirror: Simple magnetic bottle model using two parallel current loops
(coils) at z = -d and z = +d

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/plasmaMorror.py

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, tan, sin 
import time

tStart = time.time()
plt.close('all')


#%% INPUTS
q = 1.602e-19       # Particle charge (Proton, C)
m = 1.673e-27       # Particle mass (Proton, kg)
Rcoil = 0.5         # Radius of the mirror coils (m)
Zcoil = 1.0           # Distance of coils from origin (m), coils at z = -1.0 and z = +1.0
Icoil = 1.0e5       # Coil current (Amperes)
mu0 = 4*pi*1e-7
tS = 150e-6          # Simulation time (s)       

# Set Initial Conditions: Start near the center with enough perpendicular velocity
# to guide gyro-motion, and moderate parallel velocity so it bounces
# rather than escaping the loss cone.
X0 = 0.00; Y0 = 0.00; Z0 = 0  
VX0 = 0.33e5; VY0 = 0.33e5; VZ0 = 1e5                # x, y, z (m)
init_pos = [X0, Y0, Z0]
init_vel = [VX0,VY0, VZ0]             # vx, vy, vz (m/s)

#%%
# Computes a first-order approximation of an axisymmetric magnetic mirror field.
# Uses the analytical on-axis field Bz(z) and applies a Taylor expansion for Br
def get_B_field(x, y, z):
    def B_axial(z_pos):
        # Field on axis from a single loop at z_pos
        denom = (Rcoil**2 + (z - z_pos)**2)**(1.5)
        Bz = (mu0 * Icoil * Rcoil**2) / (2 * denom)
        # Derivative dBz/dz for radial field approximation
        dBz_dz = -1.5 * (mu0 * Icoil * Rcoil**2 * (z - z_pos)) / (Rcoil**2 + (z - z_pos)**2)**(2.5)
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
# Returns [vx, vy, vz, ax, ay, az] for solve_ivp
def EM(t, state):
    x, y, z, vx, vy, vz = state
    v = np.array([vx, vy, vz])
    # Get local magnetic field
    B = get_B_field(x, y, z)
    
    # Acceleration: a = (q/m) * (v x B)
    a = (q / m) * np.cross(v, B)
    
    return [vx, vy, vz, a[0], a[1], a[2]]


initial_state = init_pos + init_vel

t_span = (0, tS)                           # Simulate for 80 microseconds
t_eval = np.linspace(t_span[0], t_span[1], 4000) # 4,000 time steps for a smooth orbit

# Solve ODE
sol = solve_ivp(EM, t_span, initial_state, t_eval=t_eval, rtol=1e-8, atol=1e-10)

# Extract coordinates
t = t_eval
x, y, z = sol.y[0], sol.y[1], sol.y[2]
vx, vy, vz = sol.y[3], sol.y[4], sol.y[5]
v = (vx**2+vy**2+vz**2)**0.5
vR = (vx**2+vy**2)**0.5
xy = (x**2+y**2)**0.5


#%%  GRAPHICS  1
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (3.5,2.5))
ax.set_xlabel('x  [ m ]')
ax.set_ylabel('y  [ m ]')
ax.grid()
ax.plot(x,y,'b',lw = 2)
ax.plot(x[-1],y[-1],'ro',ms = 8)
ax.plot(x[0],y[0],'go',ms = 8)
ax.axis("equal")
fig1.tight_layout()

#%% 2
plt.rcParams['font.size'] = 12
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (3.5,2.5))
ax.set_xlabel(' t  [ $\mu$s ]')
ax.set_ylabel('x, y, R$_{xy}$  [ mm ]')
ax.grid()
ax.plot(t*1e6,x*1e3,'b',lw = 2,label = 'x')
ax.plot(t*1e6,y*1e3,'r',lw = 2,label = 'y')
ax.plot(t*1e6,xy*1e3,'g',lw = 2, label = 'R$_{xy}$')
ax.legend()
fig2.tight_layout()

#%% 3
plt.rcParams['font.size'] = 12
fig3, ax = plt.subplots(nrows=1, ncols=1,figsize = (3.5,2.5))
ax.set_xlabel(' t  [ $\mu$s ]')
ax.set_ylabel('z  [ m ]')
ax.grid()
ax.plot(t*1e6,z,'b',lw = 2)
fig3.tight_layout()

#%% 4
plt.rcParams['font.size'] = 12
fig4 = plt.figure(figsize=(3.5,2.5))
# 3D trajectory
ax = fig4.add_subplot(111, projection="3d")
f = 1e3
xP = f*x; yP = f*y; zP = f*z
ax.plot(xP, yP, zP, color="blue")
ax.plot(xP[0],yP[0],zP[0],'go',ms = 6)
ax.plot(xP[-1],yP[-1],zP[-1],'ro',ms = 6)         
ax.set_xlabel("x [mm]")
ax.set_ylabel("y [mm]")
ax.set_zlabel("z [mm]")
fig4.tight_layout()

#%%
plt.rcParams['font.size'] = 12
fig5, ax = plt.subplots(nrows=1, ncols=1,figsize = (3.5,2.5))
ax.set_xlabel(' t  [ $\mu$s ]')
ax.set_ylabel('v  [ m/s ]')
ax.grid()
ax.plot(t*1e6,vz,lw = 2,label = 'vz')
ax.plot(t*1e6,v,lw = 2,label = 'v')
ax.plot(t*1e6,vR,lw = 2,label = 'vR')
ax.legend()
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
fig5.tight_layout()


#%%
Z = linspace(-1,1,500)
dZ = Z[2] - z[1]
VR0 = sqrt(VX0**2 + VY0**2)
BZx, BZy, BZz = get_B_field(0,0,Z)
Bmax = max(BZz) 
Bmin = min(BZz)
thetaM = np.arcsin(sqrt(Bmin/Bmax))*180/pi
thetaV = np.arctan(1.414*VR0/VZ0)*180/pi
thetaV = np.arctan(VR0/VZ0)*180/pi


plt.rcParams['font.size'] = 12
fig6, ax = plt.subplots(nrows=1, ncols=1,figsize = (4,3))
ax.set_xlabel(' z  [ m ]')
ax.set_ylabel('B$_z$  [ T ]')
ax.grid()
ax.plot(Z,BZz,'b',lw = 2)
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
fig6.tight_layout()

fig7, ax = plt.subplots(nrows=1, ncols=1,figsize = (4,3))
ax.set_xlabel(' z  [ m ]')
ax.set_ylabel('B$_R$  [ $\mu$T ]')
#Q = R*1e3; ax.set_title('radius R = %0.1f mm' % Q)
ax.grid()

R = 10; BR = -0.5 * R*1e-3 * np.gradient(BZz,dZ)
ax.plot(Z,BR*1e6,lw = 2,label = R)
R = 20; BR = -0.5 * R*1e-3 * np.gradient(BZz,dZ)
ax.plot(Z,BR*1e6,lw = 2,label = R)
R = 30; BR = -0.5 * R*1e-3 * np.gradient(BZz,dZ)
ax.plot(Z,BR*1e6,lw = 2,label = R)
R = 40; BR = -0.5 * R*1e-3 * np.gradient(BZz,dZ)
ax.plot(Z,BR*1e6,lw = 2,label = R)

ax.legend(fontsize = 10)
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
fig7.tight_layout()

#%% 8
plt.rcParams['font.size'] = 10
fig8, ax = plt.subplots(nrows=1, ncols=1,figsize = (3.5,2.5))
ax.set_xlabel(' z  [ m ]')
ax.set_ylabel('v  [ m/s ]')
ax.grid()
ax.plot(z,vz,lw = 2,label = 'vz')
ax.plot(z,v,lw = 2,label = 'v')
ax.plot(z,vR,lw = 2,label = 'vR')
ax.plot(0,VR0,'go',ms = 6)
ax.plot(0,VZ0,'bo',ms = 6)
ax.legend(fontsize = 8)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
fig8.tight_layout()

#%% 9
plt.rcParams['font.size'] = 12
fig9, ax = plt.subplots(nrows=1, ncols=1,figsize = (3.5,2.5))
ax.set_xlabel(' z  [ n ]')
ax.set_ylabel('R$_{xy}$  [ mm ]')
ax.grid()
ax.plot(z,xy*1e3,'b')

fig9.tight_layout()

#%% CONSOLE DISPLAY
print('Loss cone angle = %0.2f deg' % thetaM)
print('Initial velocity: vx0 = %0.2e m/s' % VX0 +
      '  vy0 = %0.2e m/s' % VY0 + '  vz0 = %0.2e m/s' % VZ0)
print('Initial radial velocity = %0.2e' % VR0)
print('Launch angle = %0.2f deg' % thetaV)
Q = max(z); print('max z = %0.2f m' % Q)

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)

#%%
#'''
fig1.savefig('a1.png')
fig2.savefig('a2.png') 
fig3.savefig('a3.png') 
fig4.savefig('a4.png')
fig5.savefig('a5.png') 
fig6.savefig('a6.png') 
fig7.savefig('a7.png') 
fig8.savefig('a8.png') 
fig9.savefig('a9.png')  
#'''
