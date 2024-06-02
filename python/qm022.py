# -*- coding: utf-8 -*-
"""
# qm021.py        April 2024
# Ian Cooper         matlabvisualphysics@gmail.com
# QUANTUM MECHANICS
#   WAVE PARTICLE DUALITY: double slit diffraction
#  Website: https://d-arora.github.io/Doing-Physics-With-Matlab/
#  Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm021.pdf

"""

#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag
from scipy.integrate import odeint, quad, dblquad, simps

import random


tStart = time.time()


#%% Functions
def colour(wL):
    thiscolor = [0,0,0]
   # wL    = wL*1e+9    # Convert to nm.

    if wL<380: 
        thiscolor = [1,0,1]

    if (wL>=380)&(wL<440):
        thiscolor = [(440-wL)/(440-380),0,1]

    if (wL>=440)&(wL<490):
        thiscolor = [0,(wL-440)/(490-440),1]

    if (wL>=490)&(wL<510):
        thiscolor = [0,1,(510-wL)/(510-490)]

    if (wL>=510)&(wL<580):
        thiscolor = [(wL-510)/(580-510),1,0]

    if (wL>=580)&(wL<645):
        thiscolor = [1,(645-wL)/(645-580),0]

    if (wL>=645):
        thiscolor = [1,0,0]

#  The intensities fall off near limits of vision

    if wL>700:
       thiscolor = np.array(thiscolor) * (0.3 + 0.7*(780-wL)/(780-700))

    if wL<420:
       thiscolor = np.array(thiscolor) * (0.3 + 0.7*(wL-380)/(420-380))
    
    if thiscolor[0] < 0:
       thiscolor[0] = 0
    if thiscolor[1] < 0:
       thiscolor[1] = 0
    if thiscolor[2] < 0:
       thiscolor[2] = 0   
    return thiscolor

#%%
# Constants and Variables
wL = 400         # wavelength [ 400 - 700 nm]
b = 2e-4         # slit width  [ m ]
a = 3e-4         # slit separation [ m ]
L = 10e-3        # screen width  [ m ]
D = 1.0          # aperture to observation screen  [ m ]
N = 599          # grid points

#%% Computations
x = linspace(-L,L,N)       # sreen position  [ m ]
theta = x/D                # angular screen position  [rad]
k = 2*pi / (wL*1e-9)       # propagation constant  [1/m]
beta = 0.5*k*b*sin(theta)+1e-16
alpha = 0.5*k*a*sin(theta)
col = colour(wL)

# Envelope: single slit diffraction - irradiance
Is = ( sin(beta)/beta )**2

# Double slit diffraction - irradiance
Id = Is * cos(alpha)**2


#%%  GRAPHICS

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)

fig1, axes = plt.subplots(nrows=2, ncols=1)
R = 0
axes[R].set_xlabel('screen position  [ mm ]',fontsize = 12)
axes[R].set_ylabel('irradiance [ a.u. ]',color = 'black',fontsize = 12)
axes[R].xaxis.grid()
axes[R].yaxis.grid()

axes[R].plot(x*1e3,Is,lw = 0.5, color = col)

axes[R].fill_between(x*1e3, Id,color = col)

axes[R].set_xlim([-10,10])
bp = b*1e3; ap = a*1e3
axes[R].set_title('$\lambda$ = %2.0f nm  ' % wL  + 
   'b = %2.2f mm   '  % bp + 'a = %2.2f mm'  % ap , fontsize = 12, color = 'black') 

fig1.tight_layout()

R = 1
xSc, ySc = np.meshgrid(x*1e3,[0,+1]);
Is = Id**0.5
zSc = [Is,Is];
axes[R].pcolor(xSc, ySc, zSc, cmap = 'gray')
axes[R].plot(x*1e3,Is,lw = 3, color = col)
axes[R].set_axis_off()

fig1.savefig('a1.png')


#%%

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,1)
fig, ax = plt.subplots(1)

ax.set_xlim([-10,10])
ax.set_axis_off()

random.seed()
num = 500
for c in range(num):
    xP = -10+20*random.random()        # [mm]
    xQ = xP*1e-3                       # [m]
    yP = random.random()
    Q = xQ/D
    betaQ = 0.5*k*b*sin(Q)+1e-16
    alphaQ = 0.5*k*a*sin(Q)
    Is = ( sin(betaQ)/betaQ )**2
    Id = Is * cos(alphaQ)**2
    q = random.random()
    if q < Id: 
       ax.plot(xP,yP,'bo', ms = 1)
       ax.set_title('electrons = %2.0f  ' % num, fontsize = 12, color = 'black')

fig.savefig('a2.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)