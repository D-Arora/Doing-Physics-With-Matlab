# -*- coding: utf-8 -*-
"""
Oct 2026         sim001.py

VIBRATIONS

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/sim001.py

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, odeint, simpson
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, cos, sin 
import time

tStart = time.time()
plt.close('all')

#%% Define parameters
# Initial conditions: [position(0)  m, velocity(0)  m/s]
u0 = [1, 0.0]
# Time steps
N = 999
# Simulation time
t1 = 0; t2 = 12
tSpan = (t1, t2)
# Mass [kg]
m = 0.50666 
# Spring constant [N/m] 
k = 20.0  # Spring constant
# Damping coefficient [kg/s]
b = 0
# Strength of applied impulsive force
S = 10

#%% Equation of motion (EM)
def EM(t, state):
  x, v = state
  FD = 0
  if t > 3: FD = S
  if t > 6: FD = 0
  dxdt = v
  dvdt = - (k/m)*x - (b/m)*v + FD
  return [dxdt, dvdt]

#%% Solve Equation of motion using scipy.integrate.solve_ivp
t = np.linspace(tSpan[0], tSpan[1], N)
dt = t[2] - t[1]
sol = solve_ivp(EM, tSpan, u0, t_eval=t, method='RK45')

# Solution
x = sol.y[0]                   # displaement [m]
v = sol.y[1]                   # velocity    [m/s]
a = np.gradient(v,dt)          # acceleration [m/s^2]
K = 0.5*m*v**2                 # KE [J]

FR = -k*x                      # restoring force
FB = -b*v                      # damping force
Fnet = m*a                     # net force

FI = zeros(N)
FI[t>3] = S
FI[t>6] = 0


index1 = np.where(t>3)[0][0]
index2 = np.where(t>6)[0][0]-1
xI = x[index1:index2]
xImax = max(xI)
xImin = min(xI)
xE = 0.5*(xImax + xImin)

U = zeros(N)
for c in range(N):
    U[c] = 0.5*k*x[c]**2
    if t[c] > 3 and t[c] < 6: U[c]=0.5*k*(x[c]-xE)**2

E = K + U                      # Total energy [J]

#%% 1 Displacement
xP = t; yP = x
txtX = 't [ s ]'; txtY = 'x [ m ]'
plt.rcParams['font.size'] = 10
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
ax.plot(xP,yP,'b',lw = 2)
fig1.tight_layout()

#%% 2 Velocity
xP = t; yP = v
txtX = 't [ s ]'; txtY = 'v [ m.s$^{-1}$ ]'
plt.rcParams['font.size'] = 10
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
ax.plot(xP,yP,'b',lw = 2)
fig2.tight_layout()

#%% 3 Acceleration
xP = t; yP = a
txtX = 't [ s ]'; txtY = 'a [ m.s$^{-2}$ ]'
plt.rcParams['font.size'] = 10
fig3, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
ax.plot(xP,yP,'b',lw = 2)
fig3.tight_layout()

#%% 4 Phase plot
xP = x; yP = v
txtX = 'x [ m ]'; txtY = 'v [ m.s$^{-1}$ ]'
plt.rcParams['font.size'] = 10
fig4, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
ax.plot(xP,yP,'b',lw = 2)
ax.plot(xP[index1:index2],yP[index1:index2],'r',lw = 2)
ax.plot(xP[-1],yP[-1],'ro',ms=7)
ax.plot(xP[0],yP[0],'go',ms=7)


fig4.tight_layout()

#%% 5 Forces
xP = t
txtX = 't [ s ]'; txtY = 'F  [N]'
plt.rcParams['font.size'] = 10
fig5, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()

yP = Fnet;ax.plot(xP,yP,'k',lw = 1,label = 'F$_{net}$')
yP = FR; ax.plot(xP,yP,'r',lw = 1,label = 'F$_R$')
yP = FI; ax.plot(xP,yP,'b',lw = 2,label = 'F$_D$')

ax.legend()
fig5.tight_layout()

#%% 6 Energy
xP = t; yP = K
txtX = 't [ s ]'; txtY = 'K, U, E  [J]'
plt.rcParams['font.size'] = 10
fig6, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
ax.plot(xP,yP,'b',lw = 2)
yP = U; ax.plot(xP,yP,'r',lw = 1)
yP = E; ax.plot(xP,yP,'k',lw = 1)
fig6.tight_layout()

#%%
#'''
fig1.savefig('a1.png')
fig2.savefig('a2.png') 
fig3.savefig('a3.png') 
fig4.savefig('a4.png')
fig5.savefig('a5.png') 
fig6.savefig('a6.png') 
 
#'''
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)

