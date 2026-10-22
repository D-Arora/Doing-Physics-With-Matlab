# -*- coding: utf-8 -*-
"""
Oct 2026         plasma006.py

Motion of charged particles in magnetic and electric field
  [3D] motion of a proton in a constant electric field


Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/plasma001.py
"""

import numpy as np
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, tan, sin 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import time

tStart = time.time()
plt.close('all')

#%% Particle properties 
Ey = 1e6        # [v/m]      
t2 = 0.12e-6

q = 1.602e-19       # charge (C)
m = 1.67e-27         # mass (kg)
#m = 3.34e-27
# Initial position and velocity
r0 = np.array([-0.5, 0, 0.0])       # m
v0 = np.array([8e6,0,0])       # m/s

# Electric and magnetic fields
E = np.array([0.0, Ey, 0.0])        # V/m
B = np.array([0.0, 0.0, 0.0])        # T

#%%  Equations of motion
def equations(t, S):
    r = S[:3]       # position
    v = S[3:]       # velocity
    Es = E #np.array([0,0,2])
    if r[0] < -0.3: Es = np.array([0,0,0])
    if r[0] > 0.3: Es = np.array([0,0,0])
    if r[1] < -0.3: Es = np.array([0,0,0])
    if r[1] > 0.3: Es = np.array([0,0,0])
    # Lorentz force
    F = q * (Es + np.cross(v, B))
    # Acceleration
    a = F / m
    return np.concatenate((v, a))

# Initial state: [x, y, z, vx, vy, vz]
S0 = np.concatenate((r0, v0))

# Time interval
t1 = 0
N = 10000

t_eval = np.linspace(t1, t2, N)
dt = t_eval[2] - t_eval[1]

# Solve ODE
solution = solve_ivp(
    equations,
    (t1, t2),
    S0,
    t_eval=t_eval,
    rtol=1e-8,
    atol=1e-10)

# Extract position, velocities and kinetic energy
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]

vx = solution.y[3]
vy = solution.y[4]
vz = solution.y[5]

v_xy = (vx**2+vy**2)**0.5

# Acceleration in E field 
aY = q*Ey/m

ay = np.gradient(vy,dt)
aX = np.gradient(vx,dt)



#%% COSOLE SUMMARY
print("INPUTS")
print('   m = %0.2e kg' % m)
print('   q = %0.2e kg' % q)
print('   x0 = %0.2f m' % r0[0],'   y0 = %0.2f m' % r0[1], '   z0 = %0.2f m' % r0[2] )
print('   vx0 = %0.1e m/s' % v0[0],'  vy0 = %0.1e m/s' % v0[1], '  vz0 = %0.1e m/s' % v0[2])
print('   Ex = %0.2f N/C' % E[0],'   Ey = %0.2e N/C' % E[1], '   Ez = %0.2f N/C' % E[2] )
print('   Bx = %0.2f T' % B[0],'   By = %0.2f T' % B[1], '   Bz = %0.2f T' % B[2] )
print('   Simulation time t2 = %0.2e s' % t2)
print('   Time steps N = %0.0f' % N)
print('OUTPUTS')
print('   Ey ~= 0: ay = %0.1e m/s^2' % aY)

#xxx

#%%  GRAPHICS
plt.rcParams['font.size'] = 10

fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (4.3,4))
ax.set_xlabel(' x  [ m ]')
ax.set_ylabel('y  [ m ]')
ax.grid()
ax.plot([-0.5,0.5],[0,0],'m',lw =2)
ax.plot(x[0],y[0],'go',ms = 6)

ax.plot(x,y,'b',lw = 2)

ax.axis("equal")
rect = patches.Rectangle((-0.3,-0.3), 0.6,0.6, 
                         facecolor='blue', 
                         alpha=0.3,       # Transparency (0 to 1)
                         edgecolor='none') # Remove border if desired
ax.add_patch(rect)
fig1.tight_layout()

#%%
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' t  [ $\mu$s ]')
ax.set_ylabel('x,y  [ m ]')
ax.grid()
xP = t_eval*1e6; yP = x
ax.plot(xP,yP,'b',lw = 2,label = 'x')
yP = y; ax.plot(xP,yP,'r',lw = 2,label = 'y')
ax.legend()
fig2.tight_layout()

#%%
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' t  [ $\mu$s ]')
ax.set_ylabel('v$_x$, v$_y$  [ m/s ]')
ax.grid()
xP = t_eval*1e6; yP = vx
ax.plot(xP,yP,'b',lw = 2,label = 'v$_x$')
yP = vy; ax.plot(xP,yP,'r',lw = 2,label = 'v$_y$')
ax.legend()
fig3.tight_layout()

#%%
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' t  [ $\mu$s ]')
ax.set_ylabel('a$_x$, a$_y$  [ m/s$^2$ ]')
ax.grid()
xP = t_eval*1e6; yP = aX
ax.plot(xP,yP,'b',lw = 2,label = 'a$_x$')
yP = ay; ax.plot(xP,yP,'r',lw = 2,label = 'a$_y$')
ax.legend()
fig4.tight_layout()


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
'''
