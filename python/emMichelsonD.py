# -*- coding: utf-8 -*-
"""
emMichelsonD.py          mar 2025

COMPUTATIONAL OPTICS
MICHELSON INTERFEROMETER
    Planes waves: Mirror M2 tilted
    
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/emMichelson.pdf

"""

#%%
import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array 
from scipy import pi, sqrt
import matplotlib.pyplot as plt

#%%
flag = 3
if flag == 1:  wL = 700e-9; col = [1,0,0]; Col = 'Reds'     # red wavelength [m]
if flag == 2:  wL = 550e-9; col = [0,1,0]; Col = 'Greens'   # green wavelength [m]
if flag == 3:  wL = 480e-9; col = [0,0,1]; Col = 'Blues'    # blue wavelength [m]

num = 299        # Grid points
k = 2*pi/wL      # Propagtion constant  [1/m]

# Tilt angle [deg] >>>>>
Adeg = 10

A = Adeg*pi/180     # radians
dz = 1295e-9        # distance between the two mirrors

# Detector screen   [2D]
XD = 5e-6
xD = linspace(-XD,XD,num) 
zD = 1
yD = xD
xx,yy = np.meshgrid(xD,yD)


#%% Electric field and intensity at detector screen   [1D]
T =k*(zD + 2*dz) + pi
cosT = np.cos(T); sinT = sin(T)

E1 = exp(1j*T)
E2 = exp(1j*k*(zD*cos(A) + xD*sin(A))) 
E = E1 + E2
S = np.real(np.conj(E)*E)

RE1 = np.real(E1)
IE1 = np.imag(E1)


#%% Electric field and intensity at detector screen   [2D]
E1R = zeros([num,num]); E1I = zeros([num,num])
E2R = zeros([num,num]); E2I = zeros([num,num])

for c in range(num):
     E1R[:,c] = cosT
     E1I[:,c] = sinT
    
     E2R[:,c] = cos(k*(zD*cos(A) + xD[c]*sin(A))) 
     E2I[:,c] = sin(k*(zD*cos(A) + xD[c]*sin(A))) 

ER = E1R + E2R; EI = E1I +E2I 
SS = (ER**2 + EI**2)
SS = SS/amax(SS)

# Fringe separation
dxD = wL/sin(A)


#%% GRAPHICS

plt.rcParams["figure.figsize"] = (6,3.5)
fig1, ax = plt.subplots(nrows=2, ncols=1)
ax[0].set_xlabel('x$_D$  [ $\mu$ m]',fontsize = 12)
ax[0].set_ylabel('$S_D$  [ a.u ] ',fontsize = 12)
q = np.array([-5,5]); ax[0].set_xlim(q)
q = np.array([-0.2,1.2]); ax[0].set_ylim(q)

q = wL*1e9; ax[0].text(-4.5,0.1,'  $\lambda$ = %0.0f  nm' %q,fontsize = 12, 
        fontdict=dict(fontsize=15, fontweight='bold'), bbox=dict(facecolor='white',
        edgecolor='white'))
q = Adeg; ax[0].text(2.5,0.1,'  A = %0.0f  deg' %q,fontsize = 12, 
        fontdict=dict(fontsize=15, fontweight='bold'), bbox=dict(facecolor='white',
        edgecolor='white'))
q = dxD*1e6; ax[0].text(-1.2,0.90,'  $\Delta x_D$ = %0.2f  $\mu$m' %q,fontsize = 12, 
        fontdict=dict(fontsize=15, fontweight='bold'), bbox=dict(facecolor='white',
        edgecolor='white'))

ax[0].grid()
xP = xD*1e6; yP = S/max(S)
ax[0].plot(xP,yP,color = col,lw = 2)

ax[1].set_xlabel('x$_D$  [ $\mu$ m]',fontsize = 12)
ax[1].pcolor(xD*1e6,yD*1e6,SS,cmap = Col)
q = np.array([-5,5]); ax[0].set_xlim(q)
plt.yticks([])
fig1.tight_layout()


#%%
fig1.savefig('a1.png')

