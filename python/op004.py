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
              DETECTOR  YZ plane
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
N = 96
wL = 500e-9*ones(N)
A = 1e-5*ones(N)
phi = zeros(N); phi[0] = 0; phi[1] = 0
L = wL[0]

#phi[1::2] = pi
# Model parameters
c = 3e8        # speed of light in vaccum  [m/s]
k = 2*pi/wL    # propagation constant  [rad/m] 
f = c/wL       # frequency  [Hz]  
w = 2*pi*f     # angular frequency  [rad/s]
T = 1/f        # period  [s]

#%% SOURCES
s = L
xS = 0
yS = 0
zS = zeros(N)
zS[0] = -(N-1)*s/2
for c in range(N-1):
    zS[c+1] = zS[c]+s



# Detector >>>
nP = 599
Lmin = -250*L; Lmax = 250*L
x = 1000*L 
y = 0
z = linspace(Lmin,Lmax,nP)
R = zeros(nP)
U = zeros(nP)
for m in range(N):
    R = sqrt((x-xS)**2 + (y-yS)**2 + (z-zS[m])**2)
    u = A[m]*exp(1j*(k[m]*R+phi[m]))/R 
    U = U + u

I = abs(U)**2
plt.plot(z,I)
xxx


#%% Console output
wLs = 1e9*wL; Ts = 1e15*T; phis = phi/pi
print('wL [nm]   A     f [Hz]     T [fs]     phi/pi')
for m in range(2):
    print(' %2.0f ' %wLs[m] + '    %2.1e' %A[m] + '   %2.3e' %f[m]  
          + '   %  2.3f' %Ts[m] +'     %2.2f' %phis[m] ) 
print('  ')
xx = x/L; print('YZ plane  xP/wL1 = %0.0f  m' %xx)

#%% GRAPHICS
# Fig 1: wavefunction u
plt.rcParams["figure.figsize"] = (5,4)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('y  [m]',fontsize = 12)
ax.set_ylabel('z  [m]',fontsize = 12)

xP = Y; yP = Z; zP = I/np.amax(I)
IM = plt.pcolor(xP,yP,zP)

hot(); colorbar()
ax.set_yticks([Lmin,Lmin/2,0,Lmax/2,Lmax])
ax.set_aspect('equal', 'box')
ax.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
fig1.tight_layout()

  
#%% Save figures    
# fig1.savefig('a1.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)