# -*- coding: utf-8 -*-
"""
emRSBessel02.py               April 2025
 

COMPUTATIONAL OPTICS
   BESSEL BEAM
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   IZ irradiance
   x,y,z   [1D] variables   /   X,Y [2D] variables 
   

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
 Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRS02.pdf

https://dcc.ligo.org/public/0091/G1200548/001/G1200548-v1.pdf

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
#from mpl_toolkits.mplot3d import axes3d
import sympy as sym
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import jv
tStart = time.time()

#%% INPUT PARAMETERS 
# [1d]  x,y,z   [2D]  X,Y,Z
# Grid points: Q aperture space    nQ  odd number
#              P observation space nP  odd number 
NQ = 50
NP = 50

nQ = 2*NQ + 1; nP = 2*NP+1 

# Wavelength [m]
wL = 632.8e-9 

cL = 2.99792458e8     # speed of light
eps0 = 8.854187e-12   # % permittivity of free space       

#%% Aperture space: radius a / XY dimensions of aperture [m]
a = 3.182e-3
zQ = 0
xQ = linspace(-a,a,nQ)    
yQ = linspace(-a,a,nQ)
XQ, YQ = np.meshgrid(xQ,yQ)
RQ = (XQ**2 + YQ**2)**0.5

q = 16
L = a/q
EQ = cos(2*pi*RQ/L)
EQ[EQ<0] = 0
EQ[EQ>0] = 1
EQ[RQ > a] = 0
IQ = (cL*eps0/2)*np.real(np.conj(EQ)*EQ)
IQx = IQ[NQ,:]
   
#%% Observation space 
z1 = 0.5; z2 = 5
xP = 0; yP = 0;
zP = linspace(z1,z2,nP)
EP = np.zeros(nP) + np.zeros(nP)*1j      
        

#%% SETUP 
k = 2*pi/wL            # propagation constant
ik = k*1j              # jk
# Simpson [2D] coefficients
hx = 2*a/(nQ-1); hy = 2*a/(nQ-1); h = (hx*hy/9)/(2*pi)
S = np.ones(nQ)
R = np.arange(1,nQ,2);   S[R] = 4;
R = np.arange(2,nQ-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = h*scx*scy


#%% COMPUATION OF DIFFRACTION INTEGRAL FOR ELECTRIC FIELD % IRRADIANCE
for c in range(nP):
    rPQ = np.sqrt((xP - XQ)**2 + (yP - YQ)**2 + (zP[c]-zQ)**2)
    rPQ3 = rPQ*rPQ*rPQ
    kk = ik * rPQ
    MP1 = exp(kk)
    MP1 = MP1 / rPQ3
    MP2 = zP[c] * (ik * rPQ - 1)
    MP = MP1 * MP2
    EP[c] = sum(sum(EQ*MP*S))

# Irradiance +Z axis [a.u.]
cL = 2.99792458e8     # speed of light
eps0 = 8.854187e-12   # % permittivity of free space

Iz = (cL*eps0/2)*np.real(EP*np.conj(EP))

zPeak = zP[Iz==max(Iz)][0]


#%%  GRAPHICS
# 1  Aperture [1D] radial irradiance
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('$r_Q$ / a ', fontsize=12)
ax.set_ylabel('I$_Q$  ', fontsize=12)
ax.grid()
#ax.set_xlim([-1,1]); ax.set_ylim([0,1.2])
ax.plot(xQ/a,IQx,'k',lw = 2)
fig1.tight_layout()



#%%  GRAPHICS       1   z vs Iz



plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)

fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(zP, Iz,'k',lw = 2)
ax.set_xlabel('$z_P$  [m] ', fontsize=12)
ax.set_ylabel('I$_z$  [a.u.]', fontsize=12)
ax.grid()
fig2.tight_layout()
plt.show()






#%% Console summary
print('  ')
#print('Order Bessel function  n = %0.0f' %n)
print('nQ = %0.0f   '%nQ + 'nP = %0.0f' %nP)
q = wL*1e9; print('wavelength  wL = %0.0f  nm' %q )
q = a*1e3; print('aperture radius  a = %0.3f  mm' %q ) 
print('zPeak = %0.3f m   ' % zPeak) 

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


#%%
# fig1.savefig('a1.png')

 