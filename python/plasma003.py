# -*- coding: utf-8 -*-
"""
Oct 2026         plasma03.py

Motion of charged particles in magnetic and electric field
  [3D] motion of an electron in a constant magnetic field


Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/plasma001.py
"""

import numpy as np
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, tan, sin 
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import time

tStart = time.time()
plt.close('all')

#%% Particle properties 
q = 1.602e-19       # charge (C)
m = 1.67e-27        # proton mass (kg)
m = 9.11e-31        # electron mass 
# Initial position and velocity
r0 = np.array([0.2, 0, 0.0])       # m
v0 = np.array([4e6,4e6,1e6])       # m/s

# Electric and magnetic fields
E = np.array([0.0, 0.0, 0.0])        # N/C
B = np.array([0.0, 0.0, 0.4])        # T

#%%  Equations of motion
def equations(t, S):
    r = S[:3]       # position
    v = S[3:]       # velocity
    Bs = B #np.array([0,0,2])
    #if r[2] > 30: Bs = np.array([0,0,0])
    # Lorentz force
    F = q * (E + np.cross(v, Bs))
    # Acceleration
    a = F / m
    return np.concatenate((v, a))

# Initial state: [x, y, z, vx, vy, vz]
S0 = np.concatenate((r0, v0))

# Time interval
t1 = 0
t2 = 1e-6/2000
N = 10000

t_eval = np.linspace(t1, t2, N)

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

v = (vx**2+vy**2+vz**2)**0.5
K = 0.5*m*v**2

v_xy = (vx**2+vy**2)**0.5

#%% Find period / frequency / radius
peaks, _ = find_peaks(y)
peak_times = t_eval[peaks]
periods = np.diff(peak_times)
Tavg = np.mean(periods)      # [s]
favg = 1/Tavg                # [z]
R = (max(x)-min(x))/2
# Theoretical values
Bz = B[2]
RC = m*max(v_xy)/(q*Bz)     
wC = q*Bz/m            
fC = wC/(2*pi)
TC = 1/fC


#%% CONSOLE DISPLAY
print("INPUTS")
print('   m = %0.2e kg' % m)
print('   q = %0.2e kg' % q)
print('   x0 = %0.2e m' % r0[0],'   y0 = %0.2e m' % r0[1], '   z0 = %0.ef m' % r0[2] )
print('   vx0 = %0.1e m/s' % v0[0],'  vy0 = %0.1e m/s' % v0[1], '  vz0 = %0.1e m/s' % v0[2])
print('   Ex = %0.2f N/C' % E[0],'   Ey = %0.2f N/C' % E[1], '   Ez = %0.2f N/C' % E[2] )
print('   Bx = %0.2f T' % B[0],'   By = %0.2f T' % B[1], '   Bz = %0.2f T' % B[2] )
print('   Simulation time t2 = %0.2e s' % t2)
print('   Time steps N = %0.0f' % N)
print('OUTPUTS')
print('   Radius: Rsimulation = %0.2e m ' %R, '  Rtheory = %0.2e m' %RC)
print('   Cyclotron period: Tsimulation = %0.2e s ' % Tavg, '  Ttheory = %0.2e s' %TC)
print('   Cyclotron frequency: fsimulation = %0.2e Hz ' % favg, '  ftheory = %0.2e Hz' %fC)
print('   KE is constant: KE = %0.2e J' % max(K) )

#xxx

#%%  GRAPHICS
plt.rcParams['font.size'] = 10
fig1 = plt.figure(figsize=(5.5,4))
# 3D trajectory
ax1 = fig1.add_subplot(111, projection="3d")
ax1.plot(x, y, z, color="blue")
ax1.plot(x[0],y[0]*1e3,z[0]*1e3,'go',ms = 6)
ax1.plot(x[-1],y[-1],z[-1],'ro',ms = 6)         
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
fig1.tight_layout()

#%%
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' x  [ m ]')
ax.set_ylabel('y  [ m ]')
ax.grid()
ax.plot(x,y,'b',lw = 2)
ax.plot(x[0],y[0],'go',ms = 6)
ax.plot(x[-1],y[-1],'ro',ms = 6)  
ax.axis("equal")
# plt.axvspan(0,75, color='blue', alpha=0.1)
# plt.axvspan(75,95, color='red', alpha=0.1)
# plt.axvspan(95,170, color='blue', alpha=0.1)
fig2.tight_layout()

#%%
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' t  [ us ]')
ax.set_ylabel('x,y,z  [ m ]')
ax.grid()
xP = t_eval*1e6; yP = x
ax.plot(xP,yP,'b',lw = 2,label = 'x')
#yP = y; ax.plot(xP,yP,'r',lw = 2,label = 'y')
#yP = z; ax.plot(xP,yP,'m',lw = 2,label = 'z')
# plt.axvspan(0,75, color='blue', alpha=0.1)
# plt.axvspan(75,95, color='red', alpha=0.1)
# plt.axvspan(95,170, color='blue', alpha=0.1)
ax.legend()
fig3.tight_layout()
#%%

fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' t  [ us ]')
ax.set_ylabel('vx, vy,vz  [ m/s ]')
ax.grid()
xP = t_eval*1e6; yP = vx
ax.plot(xP,yP,'b',lw = 2,label = 'vx')
yP = vy; ax.plot(xP,yP,'r',lw = 2,label = 'vy')
yP = vz; ax.plot(xP,yP,'m',lw = 2,label = 'vz')
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
