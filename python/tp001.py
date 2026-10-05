# -*- coding: utf-8 -*-

"""
Oct 2026         tp001.py

Thermal conduction through a slab

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/tp001.py
"""

import numpy as np
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, tan, sin 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import simpson
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.signal import find_peaks

tStart = time.time()
plt.close('all')


#%% INPUTS
# Input flagS: building brick 1 or glass window 2
flagS = 1    # <<<<<

# Building brick
if flagS == 1:   
   nx = 121      # number of spatial elements (101)
   nt = 30000    # number of time steps (50000)
   L = 150e-3     # width
   k1 = 1.33     #  thermal conductivity W/(m.K)
   
   rho1 = 1600   #  density kg/m**3
   c1 = 900      #  specific heat capacity J/(kg.K)

# Glass window
if flagS == 2:       
   nx = 121      # number of spatial elements (101)
   nt = 30000    # number of time steps (50000)
   L = 3e-3      # width
   k1 = 0.8      #  thermal conductivity W/(m.K)
   rho1 = 2600    #  density kg/m**3
   c1 = 840     #  specific heat capacity J/(kg.K)
   

#%% SETUP
# Cross-sectional area (1 m**2)
area = 1

k = k1 * ones(nx)
c = c1 * ones(nx)
rho = rho1 * ones(nx)

TA = 5   # outside temperature  degC
TB = 25  # inside temperature   degC
    
# initial temperature (outside temperature)  degC
T = TA*np.ones((nt,nx))

# intial energy flux density  W/m**2
J = zeros((nt,nx))

# BOUNDARY CONDITIONS  degC
T[:,0] = TA
T[:,-1] = TB


#%%  SETUP
# Cells
dx = L / nx 
x = dx/2+ dx*arange(0,nx,1)

# time increment / time
dt = 0.5 * dx**2 * min(rho) * min(c) / (2 * max(k))
t = dt*arange(0,nt,1)

# Constants 
K1 = (k / (dx))
K2 = dt / (rho * c * dx)


#%% Solve differential equations
for ct in range(nt-1):
   for cx in range(1,nx):
    J[ct+1,cx] = K1[cx] * ( T[ct,cx-1] - T[ct,cx] )
    

   for cx in range(1,nx-1):
    T[ct+1,cx] = T[ct,cx] + K2[cx] * ( J[ct+1,cx] - J[ct+1,cx+1] ) 


#%% CALCULATIONS
# Energy flux dQ/dt = JA
dQ_dt = J[-1,-1] * area

# Temperature gradient dT/dx
dT_dx = (T[0,0] - T[0,-1]) / (x[-1] - x[0])
    
# Mid point temperature
b = 5; m = -dT_dx
Tcentre = m*(L/2)+ b

# time constant tau
cc = 1; indexT = round(0.5*nx+0.5)
while T[cc,indexT] < 0.63 * Tcentre:
      cc = cc + 1;
tau = t[cc] 
 
#  max time for simulation
tMax = t[-1]
 
# Temperature gradient
dT_dx = (TA - TB)/L


#%% CONSOLE DISPLAY
if flagS == 1:
    txt = 'BRICK'
    print(txt)
    print('   Simulation time = %0.0f min' %max(t/60))
    q = J[-1,-1]/1000; print('   Steady-state J = %0.2f kW/m**2' % q )
if flagS == 2:
    txt = 'GLASS WINDOW' 
    print(txt)
    print('   Simulation time = %0.0f sec' %max(t))
    q = J[-1,-1]/1000; print('   Steady-state J = %0.2f kW/m**2' % q )
print('Temperature gradient dT/dx = %0.2e' %dT_dx)    
 
#%% GRAPHICS 
def graphT(z):
    z1 = round(z);
    if flagS == 1: z2 = round(z1*dt/60)
    if flagS == 2: z2 = round(z1*dt,1)
    yP = T[z1,:]
    ax.plot(xP,yP,lw = 2,label = z2)
def graphJ(z):
    z1 = round(z);
    if flagS == 1: z2 = round(z1*dt/60)
    if flagS == 2: z2 = round(z1*dt,1)
    yP = J[z1,:]/1e3
    ax.plot(xP,yP,lw = 2,label = z2)    

zz = array([2000,5000,8000,12000,nt-1])
    
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3.4)



#%%
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' x  [ mm ]')
ax.set_ylabel('T  [ $^o$C ]')
txt = 'BRICK'
if flagS == 2: txt = 'GLASS'
ax.set_title(txt)
ax.grid()
xP = x*1e3; 

for cc in zz:
    graphT(cc)
if flagS == 1: txt = 'min'
if flagS == 2: txt = 'sec'
ax.legend(fontsize = 10,frameon=False, ncol = 1,title=txt, title_fontsize='small')
fig1.tight_layout()

#%%
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' x  [ m ]')
ax.set_ylabel('J  [ $kW.m^{-2}$ ]')
txt = 'BRICK'
if flagS == 2: txt = 'GLASS'
ax.set_title(txt)
ax.grid()
for cc in zz:
    graphJ(cc)
if flagS == 1: txt = 'min'
if flagS == 2: txt = 'sec'
ax.legend(fontsize = 10,frameon=False, ncol = 1,title=txt, title_fontsize='small')
fig2.tight_layout()


#%% temperature profile 
plt.rcParams["figure.figsize"] = (5,1.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)

[xSc, ySc] = np.meshgrid(x,[-1,+1])
zSc = zeros((2,nx))
zSc[0,:] = T[-1,:]
zSc[1,:] = T[-1,:]
plt.pcolor(xSc, ySc, zSc,cmap = 'bwr')
ax.set_yticklabels([]); plt.yticks([]); plt.xticks([]) 
fig3.tight_layout()

#%% Time evolution at centre of slab
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3.4)
fig4, ax = plt.subplots(nrows=1, ncols=1)
if flagS == 1:
   ax.set_xlabel(' t  [ m ]')
   xP = t/60
if flagS == 2:
   ax.set_xlabel(' t  [ s ]')
   xP = t   
ax.set_ylabel('T  [ $^o$C ]',color = [1,0,0])
txt = 'BRICK: slab centre'
if flagS == 2: txt = 'GLASS: slab centre'
ax.set_title(txt,fontsize = 12)

indX = round(nx/2); yP = T[:,indX]
ax.plot(xP,yP,'r',lw=2)
ax.grid()

ax2 = ax.twinx()
ax2.set_ylabel('J  [ $kW.m^{-2}$ ]',color = [0,0,1])
indX = round(nx/2); yP = J[:,indX]/1e3
ax2.plot(xP,yP,'b',lw=2)
fig4.tight_layout()

'''#%%
fig1.savefig('a1.png') 
fig2.savefig('a2.png') 
fig3.savefig('a3.png') 
fig4.savefig('a4.png') 
'''

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)