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
              DETECTOR  XY plane
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag, real
import random
from pylab import imshow,show,gray,jet,colorbar,xlabel,ylabel,title,hsv,hot

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
Lmin = -200*L; Lmax = 200*L
z = 100*L 
x = linspace(Lmin,Lmax,nP)
y = linspace(Lmin,Lmax,nP)
X, Y = np.meshgrid(x,y)
R1 = sqrt((X-x1)**2 + (Y-y1)**2 + (z-z1)**2)
R2 = sqrt((X-x1)**2 + (Y-y1)**2 + (z-z2)**2)

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
print('  ')
print(' ')
xx = 2*s/L; print('source separation 2s/wL = %0.1f' %xx)
zz = z/L; print('XY plane  zP/wL = %0.0f' %zz)

#%% GRAPHICS
# Fig 1: wavefunction u
plt.rcParams["figure.figsize"] = (5,4)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x / $\lambda$',fontsize = 12)
ax.set_ylabel('y / $\lambda$',fontsize = 12)

xP = X/L; yP = Y/L; zP = I/np.amax(I)
IM = plt.pcolor(xP,yP,zP)

hot(); colorbar()
#ax.set_yticks([Lmin,Lmin/2,0,Lmax/2,Lmax])
ax.set_aspect('equal', 'box')
#ax.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
fig1.tight_layout()

  
#%% Save figures    
# fig1.savefig('a1.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)