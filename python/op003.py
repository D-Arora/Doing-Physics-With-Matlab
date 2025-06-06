# -*- coding: utf-8 -*-

"""
op003.py    oct 2024

COMPUTATIONAL OPTICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/op003.htm


INTERFERENCE: TWO POINT SOURCES xS = 0, yS = 0, zS1 = +s, zS2 = -s

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag, real
import random
tStart = time.time()


#%%
# Inputs: wavelength [m], amplitude, initial phase angle [rad]   >>>
wL = zeros(2); wL[0] = 500e-9; wL[1] = 500e-9
A = zeros(2); A[0] = 1e-5; A[1] = 1e-5
phi = zeros(2); phi[0] = 0; phi[1] = 0

# Source  >>>
L = wL[0]
s = 10*L
x1 = 0; y1 =0; z1 = s
x2 = 0; y1 = 0; z2 = -s

# Model parameters
c = 3e8        # speed of light in vaccum  [m/s]
k = 2*pi/wL    # propagation constant  [rad/m] 
f = c/wL       # frequency  [Hz]  
w = 2*pi*f     # angular frequency  [rad/s]
T = 1/f        # period  [s]

# Detector >>>
nP = 599
xMin = 0*L; xMax = 50*L
y = 0; z = 0; 
x = linspace(xMin,xMax,nP)

R1 = sqrt((x-x1)**2 + (y-y1)**2 + (z-z1)**2)
R2 = sqrt((x-x1)**2 + (y-y1)**2 + (z-z2)**2)

u1 = A[0]*exp(1j*(k[0]*R1+phi[0]))/R1 
u2 = A[1]*exp(1j*(k[1]*R2+phi[1]))/R2

u = u1 + u2

I = abs(u)**2

#%% Console output
wLs = 1e9*wL; Ts = 1e15*T; phis = phi/pi
print('wL [nm]   A     f [Hz]     T [fs]     phi/pi')
for m in range(2):
    print(' %2.0f ' %wLs[m] + '    %2.1e' %A[m] + '   %2.3e' %f[m]  
          + '   %  2.3f' %Ts[m] +'     %2.2f' %phis[m] ) 
print(' ')
xx = 2*s/L; print('source separation 2s/wL = %0.1f' %xx)

#%% GRAPHICS
# Fig 1: wavefunction u
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x  [nm]',fontsize = 12)
ax.set_ylabel('u  [a.u.]',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
#ax.set_ylim((-Amax, Amax))
xP = x*1e9; yP = real(u)
ax.plot(xP,yP,'b',lw = 2)
fig1.tight_layout()

# Fig 2: irradiance I
plt.rcParams["figure.figsize"] = (6,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x  [nm]',fontsize = 12)
ax.set_ylabel('I  [a.u.]',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
#ax.set_ylim((-Amax, Amax))
xP = x*1e9; yP = I
ax.plot(xP,yP,'b',lw = 2)
fig2.tight_layout()

  
#%% Save figures    

fig1.savefig('a1.png')
fig2.savefig('a2.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)