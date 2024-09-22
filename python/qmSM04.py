# -*- coding: utf-8 -*-
"""
qmSM02.py    Aug 2024

QUANTUM MECHANICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSM01.pdf


STATISTICAL MECHANICS: MAXWELL - BOLTZMANN DISTRIBUTION
"""

#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np
import math
from math import factorial 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag
from scipy.integrate import odeint, quad, dblquad, simps

tStart = time.time()


#%% INPUTS >>>>   
# Temperature T [K]
T = 300
# amu [kg] / molecular mass of gas
M = 32
amu = 1.6605e-27
# mass of gas molecule m [kg]  = amu * molecule mass of gas
m = amu*M
# Velocity range  [m/s]
vMin = 0; vMax = 5000
# Velocity limits for probability calculation [m/s]
v1= 000; v2 = 1000
# Grid size
num = 9999

#%% SETUP
# Boltzmann constant
kB = 1.38e-23
# Universal gas constant 
R = 8.31
v = linspace(vMin, vMax, num)


#%% Maxwell speed distribution
fM = 4*pi * (m/(2*pi*kB*T))**1.5 *  v**2 * exp(-(m*v**2) / (2*kB*T))
# CHECK: Normalization: area = 1: may need to ajust vMax and num 
area = simps(fM,v)

# Average speed
fn = v*fM
vAvg = simps(fn,v)
# Most probable speed
target_value = max(fM)
index1 = np.argmin(np.abs(fM - target_value))
vP = v[index1]
# RMS speed
fn = v**2 * fM
vRMS =sqrt(simps(fn,v))

# Theoretical values
vTavg = sqrt(8*kB*T/(pi*m))
vTp   = sqrt(2*kB*T/m)
vTrms = sqrt(3*kB*T/m)
print(vTavg)
print(vTp)
print(vTrms)


#%%  Probability of finding particle between v1 and v2
target_value = v1
a = np.argmin(np.abs(v - target_value))
target_value = v2
b = np.argmin(np.abs(v - target_value))

prob = simps(fM[a:b],v[a:b])


#%% CONSOLE OUTPUT
print('  ')
print('Numerical calculations')
print('T = %0.0f K' %T + '   M = %0.0f kg' %M)
print('Normalization area under distribution curve = 1:  area = %0.7f' %area)
print('Most probable speed: vP = %0.0f m/s' %vP)
print('average speed: vAvg = %0.0f m/s' %vAvg)
print('RMS speed: vRMS = %0.0f m/s' %vRMS)
print('Probability of finding particle v1 < v < v2')
print('v1 = %0.0f' % v1 + '   v2 = %0.0f' % v2 +    '   prob = %0.3f' % prob)


print('  ')
print('Theoretical calculations')
print('most probable: vP = %0.0f m/s' %vTp)
print('average: vAvg = %0.0f m/s' %vTavg)
print('RMS: vRMS = %0.0f m/s' %vTrms)



#%% Figure 1: Plot of Maxwell speed distribution
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('v [ m/s ]',fontsize = 12)
ax.set_ylabel('f$_M$  [ x10$^{-3}$  s/m ]',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
xP = v
yP = fM*1e3
ax.plot(xP,yP,'b',lw = 2)

xP = [vP,vP]; yP = [0,max(fM)*1e3]
ax.plot(xP,yP,'r',lw = 1, label = 'vP')

xP = [vAvg,vAvg]
target_value = vAvg
index1 = np.argmin(np.abs(v - target_value))
y = fM[index1]
yP = [0,y*1e3]
ax.plot(xP,yP,'k',lw = 1, label = 'vAvg')

xP = [vRMS,vRMS]
target_value = vRMS
index1 = np.argmin(np.abs(v - target_value))
y = fM[index1]
yP = [0,y*1e3]
ax.plot(xP,yP,'m',lw = 1, label = 'vRMS')
ax.legend()

ax.set_title('T = %0.0f K' %T + '    m = %0.2e kg' %m, fontsize = 10)
ax.set_ylim([0,2.5])

fig1.tight_layout()


#%% Figure 2: Plot of Maxwell speed distribution
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('v [ m/s ]',fontsize = 12)
ax.set_ylabel('f$_M$  [ x10$^{-3}$  s/m ]',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
xP = v
yP = fM*1e3
ax.plot(xP,yP,'b',lw = 3)

ax.fill_between(v[a:b], 1e3*fM[a:b],color = [0.8,0.2,0.2],alpha=0.2)

ax.set_title('T = %0.0f K' %T + '  m = %0.2e kg' %m +'   prob = %0.3f '%prob,
             fontsize = 10)

fig2.tight_layout()


#%% SAVE FIGURES
fig1.savefig('a1.png')
fig2.savefig('a2.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)



