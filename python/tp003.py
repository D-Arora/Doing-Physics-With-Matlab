# -*- coding: utf-8 -*-

"""
Oct 2026         tp003.py

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
nt = 500000    # number of time steps (50000)
# Building brick 1 / insulation fibre 2
# number of spatial elements 
#nx1 = 81
#nx2 = 81
nx  = 170       
# width  [m]  
L1 = 2*75e-3 
L2 = 20e-3
L = L1+ L2
#  thermal conductivity W/(m.K)     
k1 = 1.33     
k2 = 0.038
#  density kg/m**3  
rho1 = 1600  
rho2 = 200
#  specific heat capacity J/(kg.K)
c1 = 900
c2 = 2100      
# Glass window

#%% SETUP
# Cross-sectional area (1 m**2)
area = 1

k = k1 * ones(nx)
c = c1 * ones(nx)
rho = rho1 * ones(nx)
k[75:95] = k2
c[75:95]  = c2
rho[75:95] = rho2

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

# Centre temperature
Tc = ( TB+(k1/k2)*TA )/( 1 + k1/k2 )

#  max time for simulation
tMax = t[-1]
 

#%% CONSOLE DISPLAY

q = max(t)/3600; print('   Simulation time = %0.0f hrs' % q)
q = J[-1,-1]; print('   Steady-state J = %0.2f W/m**2' % q )
q = Tc; print('   Tcentre = %0.1f degC' % q )
 

#%%    
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3.4)

fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' x  [ mm ]')
ax.set_ylabel('T  [ $^o$C ]')
#ax.text(28,22,'BRICK')
#ax.text(103,22,'FIBRE')
#ax.text(43,8,'T$_C$ = %0.1f $^o$C'%Tc)
ax.set_xlim([-0,170])
ax.grid()
xP = x*1e3; yP = T[-1,:]
ax.plot(xP,yP,'b',lw = 2) 
plt.axvspan(0,75, color='blue', alpha=0.1)
plt.axvspan(75,95, color='red', alpha=0.1)
plt.axvspan(95,170, color='blue', alpha=0.1)
fig1.tight_layout()

#%%
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' x  [ m ]')
ax.set_ylabel('J  [ $kW.m^{-2}$ ]')
txt = 'BRICK'
#ax.set_title(txt)
ax.grid()
J[-1,0] = J[-1,1]
xP = x*1e3; yP = J[-1,:]
#ax.set_ylim([-12,0])
ax.plot(xP,yP,'r',lw=2)
# for cc in zz:
#     graphJ(cc)
# if flagS == 1: txt = 'min'
# if flagS == 2: txt = 'sec'
# ax.legend(fontsize = 10,frameon=False, ncol = 1,title=txt, title_fontsize='small')
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

#%% Time evolution temperature centre of brick and fibre
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3.4)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' t  [ hrs ]')
xP = t/3600
ax.set_ylabel('T  [ $^o$C ]')
yP = T[:,30]
ax.plot(xP,yP,'r',lw=2,label ='centre brick')
yP = T[:,90]
ax.plot(xP,yP,'b',lw=2,label ='centre fibre')
ax.grid()
ax.legend(fontsize = 10,frameon=False, ncol = 1)
fig4.tight_layout()

#%% Time evolution aenergy flux density centre of brick and fibre
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3.4)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(' t  [ hrs ]')
xP = t/3600
ax.set_ylabel('J  [ W.m$^{-2}$ ]')
yP = J[:,30]
ax.plot(xP,yP,'r',lw=2,label ='centre brick')
yP = J[:,90]
ax.plot(xP,yP,'b',lw=2,label ='centre fibre')
ax.grid()
ax.legend(fontsize = 10,frameon=False, ncol = 1)
fig5.tight_layout()

'''
fig1.savefig('a1.png') 
fig2.savefig('a2.png') 
fig3.savefig('a3.png') 
fig4.savefig('a4.png') 
fig5.savefig('a5.png') 
'''


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)