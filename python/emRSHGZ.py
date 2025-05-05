# -*- coding: utf-8 -*-
"""
emRSHGZ.py               30 April 2025
 
 HERMITE-GAUSS FUNCTION
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   IZ irradiance
   x,y,z   [1D] variables   /   X,Y [2D] variables 
   
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRSHG.pdf
"""

import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, sqrt 
import matplotlib.pyplot as plt
import time
from scipy import special

tStart = time.time()

plt.close('all')

#%% INPUT PARAMETERS 
# Grid points: Q aperture space   /   P observation space 
NQ = 100
NP = 100

nQ = 2*NQ + 1; nP = 2*NP+1 

# Wavelength [m]
wL = 632.8e-9  
 
#%% APERTURE SPACE
#  HERMITRE-GAUSSS [2D] FUNCTION
m = 1           # x and y orders
n = 1            
a = 4e-4        # radius of circular aperture [m]
w0 = 0.5*a      # standard eviation of Gausssian function
q = sqrt(2)/w0  # constant

xQ = np.linspace(-a, a, nQ)
yQ = xQ
XQ,YQ = np.meshgrid(xQ,yQ) 
RQ = (XQ**2 + YQ**2)**0.5   
# Gaussian function
E0 = exp(-RQ**2/(2*w0**2))
# Hermite functions
xH = special.hermite(m, monic=False)
yH = special.hermite(n, monic=False)
Hx  = xH(q*XQ); Hy = yH(q*YQ)
HH = Hx*Hy*E0
EQ = HH
IQ = np.real(np.conj(EQ)*EQ)         #   IQ(x,y,0)
IQx = IQ[NQ,:]                       #   IQ(x,y=0,0)

#%% OBSERVATION SPACE   [distances m]
xP = 0; yP = 0;
z1 = 0.01; z2 = 1; zP = linspace(z1,z2,nP)       
EP = np.zeros(nP) + np.zeros(nP)*1j    

#%% SETUP 
k = 2*pi/wL            # propagation constant
ik = k*1j              # jk
# Simpson's [2D] coefficients
S = np.ones(nQ)
R = np.arange(1,nQ,2);   S[R] = 4;
R = np.arange(2,nQ-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy

#%% COMPUATION OF DIFFRACTION INTEGRAL FOR ELECTRIC FIELD and IRRADIANCE
for c1 in range(nP):
    rPQ = np.sqrt((xP - XQ)**2 + (yP - YQ)**2 + zP[c1]**2)
    rPQ3 = rPQ*rPQ*rPQ
    kk = ik * rPQ
    MP1 = exp(kk)
    MP1 = MP1 / rPQ3
    MP2 = zP[c1] * (ik * rPQ - 1)
    MP = MP1 * MP2
    EP[c1] = sum(sum(EQ*MP*S))

# Intensity (irradiance) optical axis [a.u.]
Iz = np.real(EP*np.conj(EP))
Iz = Iz/amax(amax(Iz))

#%%  GRAPHICS
# 1   Aperture IQ(x,y=0)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(xQ/a, IQx,'k',lw = 2)

ax.set_title('m = %0.0f     '%m + 'n = %0.0f' %n,fontsize = 10)
ax.set_xlabel('$x_Q$/a ', fontsize=12)
ax.set_ylabel('IQ$_z$  [a.u.]', fontsize=12)
ax.grid()
fig1.tight_layout()

# 2   Aperture  IQ(x,y)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
# fig2.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.0,\
#                     right = 1, hspace = 0.10,wspace=0.1)   
cf = ax.pcolor(XQ,YQ,EQ**2, cmap='jet')  
    # Greys_r
ax.set_aspect('equal', adjustable='box')
ax.set_title('m = %0.0f     '%m + 'n = %0.0f' %n,fontsize = 10)
ax.set_xticks([]); ax.set_yticks([])
fig2.colorbar(cf, ax=ax)
fig2.tight_layout()

# 3   Observation space  IQ(x = 0, y = 0, z)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(zP, Iz,'k',lw = 2)

ax.set_xlim([0,z2]); ax.set_ylim([0,1.1]);
ax.set_xlabel('$z_P$  [m] ', fontsize=12)
ax.set_ylabel('I$_z$  [a.u.]', fontsize=12)
ax.grid()
fig3.tight_layout()

#%% Console summary
print('  ')
print('nQ = %0.0f   '%nQ + 'nP = %0.0f' %nP)
q = wL*1e9; print('wavelength  wL = %0.0f  nm' %q )
q = a*1e3; print('aperture radius  a = %0.3f  mm' %q ) 
print('z1 = %0.2f m   '%z1 + 'z2 = %0.2f' %z2 )  
print('m = %0.0f   '%m + 'n = %0.0f' %n) 

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

#%%

# fig1.savefig('a1.png')
# fig2.savefig('a1.png')
# fig3.savefig('a1.png')