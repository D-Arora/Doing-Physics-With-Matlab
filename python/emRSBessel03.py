# -*- coding: utf-8 -*-
"""
emRSBessel03.py               April 2025
 

COMPUTATIONAL OPTICS
   BESSEL BEAM
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   ZX irradiance
   x,y,z   [1D] variables   /   X,Y [2D] variables 
   

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
 Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRSBessel.pdf


https://www.researchgate.net/publication/232992585_Bessel_beams_Diffraction_in_a_new_light

"""

import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax 
from numpy.linalg import eig
from scipy.integrate import odeint, quad, dblquad, simps
from scipy import pi, sqrt
from scipy.special import j1
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import time
from scipy.special import jv
import sympy as sym
from mpl_toolkits.mplot3d import Axes3D

tStart = time.time()

#%% INPUT PARAMETERS 
# [1d]  x,y,z   [2D]  X,Y,Z
# Grid points: Q aperture space    nQ format odd number
#              P observation space nP format integer * 4 + 1
nQ = 199           # 199
num = 59          # 59       
nP = num*4+1
         
# Wavelength [m]
wL = 632.8e-9 

# Order of Bessel function and max r value
n = 4.5; rMax = 15
       
#%% Aperautre space
a = 1e-3        # radius  [m]
zQ = 0
xQ = linspace(-a,a,nQ)    
yQ = linspace(-a,a,nQ)
XQ, YQ = np.meshgrid(xQ,yQ)
RQ = (XQ**2 + YQ**2)**0.5
EQ = jv(n, RQ*rMax/a)
EQ[RQ > a] = 0
IQ = np.real(np.conj(EQ)*EQ)


#%% Observation space 
yP = 0     # <<<<<<
Lx = 1*a

Lz1 = 0.05
Lz2 = 1

xP = linspace(-Lx,Lx,nP)     
zP = linspace(Lz1,Lz2,nP)
ZP, XP = np.meshgrid(zP,xP)  
EP = np.zeros([nP,nP])+np.zeros([nP,nP])*1j         


#%% SETUP 
k = 2*pi/wL            # propagation constant
ik = k*1j              # jk
# Simpson [2D] coefficients
S = np.ones(nQ)
R = np.arange(1,nQ,2);   S[R] = 4;
R = np.arange(2,nQ-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy


#%% COMPUATION OF DIFFRACTION INTEGRAL FOR ELECTRIC FIELD % IRRADIANCE
for c1 in range(nP):
    for c2 in range(nP):
         rPQ = np.sqrt((XP[c1,c2] - XQ)**2 + (yP - YQ)**2 + (ZP[c1,c2] - zQ)**2)
         rPQ3 = rPQ*rPQ*rPQ
         kk = ik * rPQ
         MP1 = exp(kk)
         MP1 = MP1 / rPQ3
         MP2 = ZP[c1,c2] * (ik * rPQ - 1)
         MP = MP1 * MP2
         EP[c1,c2] = sum(sum(EQ*MP*S))

# Intensity (irradiance) ZX plane [a.u.]
IZX = np.real(EP*np.conj(EP))
IZX = IZX/amax(IZX)


#%%  GRAPHICS
# 1  Aperture intensity 
plt.rcParams["figure.figsize"] = (5,4)
fig1 = plt.figure()
ax = fig1.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([]); ax.set_xticks([]); ax.set_yticks([])
ax.set(zlabel=None)
ax.set_title('J$_{%0.1f}$' %n,fontsize = 14)
ax.plot_surface(XQ/a,YQ/a,IQ,cmap = 'jet')
ax.view_init(elev=57, azim=-57, roll=0)
fig1.tight_layout()

#%%  2 Aperture intensity 
plt.rcParams["figure.figsize"] = (4,4)
fig2, ax = plt.subplots(nrows=1, ncols=1)
plt.pcolor(XQ/a,YQ/a,IQ, cmap='Greys_r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('$x_Q$ / a ', fontsize=12)
ax.set_ylabel('$y_Q$ / a ', fontsize=12)
ax.set_title('J$_{%0.1f}$' %n,fontsize = 14)
ax.set_xticks(np.arange(-1,1.4,0.5))
ax.set_yticks(np.arange(-1,1.4,0.5))
fig2.tight_layout()

#%%   3 aperture radial intensity
xZero = int(np.floor(nQ/2))
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(xQ/a,IQ[xZero,:],'k',lw = 2)
ax.set_xlabel('$r_Q$ / a ', fontsize=12)
ax.set_ylabel('I$_Q$  [a.u.]', fontsize=12)
ax.set_title('J$_{%0.1f}$' %n,fontsize = 14)
ax.grid()
fig3.tight_layout()

#%%  4 Observation space intensity 
plt.rcParams["figure.figsize"] = (6,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
plt.pcolor(ZP,XP/a,IZX**0.4, cmap='jet')   #Greys_r
ax.set_xlabel('$z_P$ [ m ]', fontsize=12)
ax.set_ylabel('x$_P$ / a', fontsize=12)        
fig4.tight_layout()


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.0f' %tExe)


#%% Console summary
print('  ')
print('nQ = %0.0f   ' %nQ  + 'nP = %0.0f'%nP)
q = a*1e3; print('aperture radius = %0.3f  mm  ' %q)
q = wL*1e9; print('Wavelength wL = %0.1f  nm  ' %q)
print('Execution time %0.0f  s ' %tExe)


#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png') 


