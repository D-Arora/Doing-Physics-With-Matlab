# -*- coding: utf-8 -*-
"""
emMichelsonA1.py          mar 2025

COMPUTATIONAL OPTICS
MICHELSON INTERFEROMETER
   plane wave - screen intensity
 
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/emMichelson.pdf
    
"""


#%%
import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array 
import matplotlib.pyplot as plt


#%%
flag = 1
if flag == 1:  wL = 700e-9; col = [1,0,0]     # red wavelength [m]
if flag == 2:  wL = 550e-9; col = [0,1,0]    # green wavelength [m]
if flag == 3:  wL = 480e-9; col = [0,0,1]    # blue wavelength [m]

num = 999        # Grid points
k = 2*pi/wL      # Propagtion constant  [1/m]
L = 1000e-9      # z length factor  [m]
# Sources 1 and 2
x1 = 0; y1 = 0; z1 = 0
x2 = 0; y2 = 0
z2 = linspace(-L,L, num)

# Detector screen
xD = 0; yD = 0; zD = 100*L


#%% Electric field and intensity at detector screen
E1 = exp(1j*k*zD)
E2 = exp(1j*(k*(zD+2*z2)+pi))
E = E1 + E2
S = np.real(np.conj(E)*E)

#%%
plt.rcParams["figure.figsize"] = (6,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('z$_2$  [ $\mu$m]  ',fontsize = 12)
ax.set_ylabel('$S_D$  [ a.u ] ',fontsize = 12)
q = np.array([-1,1]); ax.set_xlim(q)
#q = np.arange(-2,0,1);ax.set_xticks(q)
q = wL*1e9; ax.text(-0.95,0.88,'  $\lambda$ = %0.0f  nm' %q,fontsize = 12,
        fontdict=dict(fontsize=15, fontweight='bold'), bbox=dict(facecolor='white',
        edgecolor='white'))
q = zD*1e3; ax.text(-0.95,0.7,'  $z_D$ = %0.2f  mm' %q,fontsize = 12,
        fontdict=dict(fontsize=15, fontweight='bold'), bbox=dict(facecolor='white',
        edgecolor='white'))
ax.grid()
xP = z2*1e6; yP = S/max(S)
ax.plot(xP,yP,color = col,lw = 2)
fig1.tight_layout()

#%%
fig1.savefig('a1.png')
