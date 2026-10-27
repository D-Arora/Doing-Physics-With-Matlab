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
u0 = [0, 0.0]
# Time steps
N = 999
# Simulation time
t1 = 0; t2 = 10
tSpan = (t1, t2)
# Mass [kg]
m = 0.50666 
# Spring constant [N/m] 
k = 20.0  # Spring constant
# Damping coefficient [kg/s]
b = 0
# Amplitude of the driving force [N] 
A = 2 
# Frequency of the driving force [rad/s]
fD = 1.02
wD = 2*pi*fD  

#%% Equation of motion (EM)
def EM(t, state):
  x, v = state
  dxdt = v
  dvdt = - (k/m)*x - (b/m)*v + (A/m)*cos(wD*t)
  return [dxdt, dvdt]

#%% Solve Equation of motion using scipy.integrate.solve_ivp
t = np.linspace(tSpan[0], tSpan[1], N)
sol = solve_ivp(EM, tSpan, u0, t_eval=t, method='RK45')

# Solution
x = sol.y[0]                   # displaement [m]
v = sol.y[1]                   # velocity    [m/s]
FD = (A/m)*cos(wD*t)          # driving force [N]
a = - (k/m)*x - (b/m)*v + FD   # acceleration [m/s^2]
K = 0.5*m*v**2                 # KE [J]
U = 0.5*k*x**2                 # PE [J]
E = K + U                      # Total energy [J]
FR = -k*x                      # restoring force
FB = -b*v                      # damping force
Fnet = m*a                     # net force

#%% FOURIOR ANALYSIS
p  = 1j*2*pi
Fmax = 3
Fmin = 0
nF = 999
F = np.linspace(Fmin,Fmax,nF)
HR = np.zeros(nF); HI = HR; H = HR; H = HR+HR*1j
for c in range(nF):
     g = x*np.exp(p*F[c]*t)
     H[c] = simpson(g,t)

psd = 2*np.conj(H)*H   
psd = np.real(psd/max(psd))   

Find = np.where(np.real(psd)>= 0.98)
II = Find[-1]
II = II[-1]
Fpeak = abs(F[II])

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
ax.plot(xP[0],yP[0],'go',ms=6)
fig4.tight_layout()

#%% 5 Forces
xP = t
txtX = 't [ s ]'; txtY = 'F  [N]'
plt.rcParams['font.size'] = 10
fig5, ax = plt.subplots(nrows=1, ncols=1,figsize = (3.5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()

yP = Fnet;ax.plot(xP,yP,lw = 1,label = 'Fnet')
yP = FD; ax.plot(xP,yP,lw = 1,label = 'FD')
yP = FB; ax.plot(xP,yP,lw = 1,label = 'FB')
yP = FR; ax.plot(xP,yP,lw = 1,label = 'FR')
yP = Fnet;ax.plot(xP,yP,lw = 2,label = 'Fnet')

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
yP = U; ax.plot(xP,yP,'r',lw = 2)
yP = E; ax.plot(xP,yP,'k',lw = 2)
fig6.tight_layout()

#%% 7 Spectrum
xP = F; yP = psd
txtX = 'f [ Hz ]'; txtY = 'psd'; txt = ('f$_{peak} = $ %0.2f Hz' %Fpeak)
plt.rcParams['font.size'] = 10
fig7, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.set_title(txt)
ax.grid()
ax.plot(xP,yP,'b',lw = 2)
fig7.tight_layout()

#%% MOTION GRAPHS
plt.rcParams['font.size'] = 10
fig8, ax = plt.subplots(nrows=3, ncols=2,figsize = (6,8))

R = 0; C = 0
xP = t; yP = x
txtX = 't [ s ]'; txtY = 'x  [ m ]'
ax[R,C].set_xlabel(txtX)
ax[R,C].set_ylabel(txtY)
ax[R,C].grid()
ax[R,C].plot(xP,yP,'b',lw = 2)

R = 0; C = 1
xP = t; yP = v
txtX = 't [ s ]'; txtY = 'v  [ m.s$^{-1}$ ]'
ax[R,C].set_xlabel(txtX)
ax[R,C].set_ylabel(txtY)
ax[R,C].grid()
ax[R,C].plot(xP,yP,'b',lw = 2)

R = 1; C = 0
xP = x; yP = v
txtX = 'x [ m ]'; txtY = 'v  [ m.s$^{-1}$ ]'
ax[R,C].set_xlabel(txtX)
ax[R,C].set_ylabel(txtY)
ax[R,C].grid()
ax[R,C].plot(xP,yP,'b',lw = 2)
ax[R,C].plot(xP[0],yP[0],'go',ms=6)

R = 1; C = 1
xP = t; yP = Fnet
txtX = 't [ s ]'; txtY = 'F  [N]'
plt.rcParams['font.size'] = 10
ax[R,C].set_xlabel(txtX)
ax[R,C].set_ylabel(txtY)
ax[R,C].plot(xP,yP,lw = 2,label = 'Fnet')
yP = FD; ax[R,C].plot(xP,yP,lw = 2,label = 'FD')
yP = FB; ax[R,C].plot(xP,yP,lw = 2,label = 'FB')
yP = FR; ax[R,C].plot(xP,yP,lw = 2,label = 'FR')
ax[R,C].legend()
ax[R,C].grid()


R = 2; C = 0
xP = t; yP = K
txtX = 't [ s ]'; txtY = 'energy  [ J ]'
plt.rcParams['font.size'] = 10
ax[R,C].set_xlabel(txtX)
ax[R,C].set_ylabel(txtY)

ax[R,C].plot(xP,yP,lw = 2,label = 'K')
yP = U; ax[R,C].plot(xP,yP,lw = 2,label = 'U')
yP = E; ax[R,C].plot(xP,yP,lw = 2,label = 'E')
ax[R,C].legend()
ax[R,C].grid()

R = 2; C = 1
xP = F; yP = psd
txtX = 'f [ Hz ]'; txtY = 'psd'; txt = ('f$_{peak} = $ %0.2f Hz' %Fpeak)
ax[R,C].set_title(txt)
plt.rcParams['font.size'] = 10
ax[R,C].set_xlabel(txtX)
ax[R,C].set_ylabel(txtY)
ax[R,C].plot(xP,yP,lw = 2)
ax[R,C].grid()

fig8.tight_layout()

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
#fig9.savefig('a9.png')  
#'''

