# -*- coding: utf-8 -*-
"""
emRS02.py               April 2025
 

COMPUTATIONAL OPTICS
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   IZ irradiance
   x,y,z   [1D] variables   /   X,Y [2D] variables 
   

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
 Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRS02.pdf
"""

import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, sqrt 
from numpy.linalg import eig
from scipy.integrate import odeint, quad, dblquad, simps
#from scipy import pi, sqrt
from scipy.special import j1
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import time
#from mpl_toolkits.mplot3d import axes3d
import sympy as sym
from mpl_toolkits.mplot3d import Axes3D

tStart = time.time()

#%% INPUT PARAMETERS 
# [1d]  x,y,z   [2D]  X,Y,Z
# Grid points: Q aperture space    nQ  odd number
#              P observation space nP  odd number 
NQ = 100
NP = 100

nQ = 2*NQ + 1; nP = 2*NP+1 
         
# Wavelength [m]
wL = 632.8e-9 
       
# Aperutre space: radius a / XY dimensions of aperture [m]
a = 4e-4
aQx = 2*a; aQy = 2*a           
 
# Observation spacehalf: h  [m] 
xP = 0; yP = 0;
z1 = 0.005; z2 = 0.35
zP = linspace(z1,z2,nP)       
         

#%% Aperture space
xQmin = -aQx/2;  xQmax = aQx/2
yQmin = -aQy/2;  yQmax = aQy/2
xQ = linspace(xQmin,xQmax,nQ)    
yQ = linspace(yQmin,yQmax,nQ)
XQ, YQ = np.meshgrid(xQ,yQ)

RQ = (XQ**2 + YQ**2)**0.5
EQ = np.ones([nQ,nQ])
EQ[RQ > xQmax] = 0
IQ = np.real(np.conj(EQ)*EQ)

#%% Observation space
EP = np.zeros(nP) + np.zeros(nP)*1j         

#%% Simpson [2D] coefficients
S = np.ones(nQ)
R = np.arange(1,nQ,2);   S[R] = 4;
R = np.arange(2,nQ-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy

#%% COMPUATION OF DIFFRACTION INTEGRAL FOR ELECTRIC FIELD % IRRADIANCE
k = 2*pi/wL            # propagation constant
ik = k*1j              # jk

for c1 in range(nP):
    rPQ = np.sqrt((xP - XQ)**2 + (yP - YQ)**2 + zP[c1]**2)
    rPQ3 = rPQ*rPQ*rPQ
    kk = ik * rPQ
    MP1 = exp(kk)
    MP1 = MP1 / rPQ3
    MP2 = zP[c1] * (ik * rPQ - 1)
    MP = MP1 * MP2
    EP[c1] = sum(sum(EQ*MP*S))

# Intensity (irradiance) XY plane [a.u.]
Iz = np.real(EP*np.conj(EP))
Iz = Iz/amax(amax(Iz))

# Rayleigh length
RL = 4*a**2/wL


#%% Console summary
print('  ')
print('nQ = %0.0f   '%nQ + 'nP = %0.0f' %nP)
q = wL*1e9; print('wavelength  wL = %0.0f  nm' %q )
q = a*1e3; print('aperture radius  a = %0.3f  mm' %q ) 
print('z1 = %0.2f m   '%z1 + 'z2 = %0.2f' %z2 )  
#q = Iz[-1]/Iz[0]; print('Iz(z1) = %0.3f   '%Iz[0] + 'Iz(z2) = %0.3f' %Iz[-1] +
#      '   Iz(z2)/Iz(z1) = %0.3f' %q) 
print('Rayleigh length  RL wL = %0.2f  nm' % RL )


#%%  GRAPHICS       1   z vs Iz
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)

fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(zP, Iz,'k',lw = 2)
ax.plot([RL,RL], [0,1],'r',lw = 1)
ax.set_xlabel('$z_P$  [m] ', fontsize=12)
ax.set_ylabel('I$_z$  [a.u.]', fontsize=12)
ax.set_title('RL = %0.2f  m' %RL)
ax.set_xlim([0,z2])
ax.grid()
fig1.tight_layout()


#%%
fig1.savefig('a1.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
