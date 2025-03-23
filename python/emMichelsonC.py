# -*- coding: utf-8 -*-

"""
emMichelsonC.py

COMPUTATIONAL OPTICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emMichelson.pdf


INTERFERENCE: TWO POINT SOURCES xS = 0, yS = 0, zS1 = +s, zS2 = -s
              Detector  XY plane: circular fringe pattern 
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag, real



#%%
flag = 1
if flag == 1:  wL = 700e-9; col = 'Reds'     # red wavelength [m]
if flag == 2:  wL = 550e-9; col = 'Greens'    # green wavelength [m]
if flag == 3:  wL = 480e-9; col = 'Blues'    # blue wavelength [m]

N = 199        # Grid points

k = 2*pi/wL    # Propagation constant  [1/m]

# Detector space - observation screen
L = 700e-9
zD = 60*L
XD = 80*L
xD = linspace(-XD,XD,N)
yD = linspace(-XD,XD,N)

X, Y = np.meshgrid(xD,yD)

# Two point sources: Source 1 at Origin (0,0,0) 
x1 = 0; x2 = 0
y1 = 0; y2 = 5*L
z1 = 0; z2 = -10.0*L

# Calculation of electric fields and screen intensity S  [a.u.]
r1 = ( (X - x1)**2 + (Y - y1)**2 + (zD-z1)**2 )**0.5
r2 = ( (X - x2)**2 + (Y - y2)**2 + (zD-z2)**2 )**0.5
E1 = exp(1j*k*r1)/r1
E2 = exp(1j*(k*r2+pi))/r2
E = E1+E2
S = np.real(np.conj(E)*E)


#%% GRAPHICS
plt.rcParams["figure.figsize"] = (4.2,2.8)
fig1, ax = plt.subplots(nrows=1, ncols=1)
                        
xP = X/L; yP = Y/L; zP = S/np.amax(S)
plt.pcolor(xP,yP,zP)

cmap = plt.get_cmap(col)
plt.set_cmap(cmap)
#colorbar()
ax.set_aspect('equal', 'box')
plt.xticks([])
plt.yticks([])

fig1.tight_layout()

  
#%% Save figures    
fig1.savefig('a1.png')


