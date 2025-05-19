# -*- coding: utf-8 -*-
"""
emRSFBZ.py               18 May 2025
 
 FOCUSED BEAM  focal length f
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   ANNULAR APERTURE
   IZ irradiance
   x,y,z   [1D] variables   /   X,Y [2D] variables 
   
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRSFB01.pdf
"""

import numpy as np
from numpy import pi, arange, exp, sin, linspace, zeros, amax, sqrt, real
import matplotlib.pyplot as plt
import time
from scipy import special
from scipy.signal import find_peaks
from scipy.special import jv
from mpl_toolkits.mplot3d import Axes3D
tStart = time.time()

plt.close('all')

#%% INPUT PARAMETERS 
# Grid points: Q aperture space   /   P observation space 
NQ = 99
NP = 99

nQ = 2*NQ + 1; nP = 2*NP+1 

sf = 0.25          # scaling factor for IXY plots

# Wavelength [m]
wL = 500e-9  

k = 2*pi/wL            # propagation constant
ik = k*1j              # jk
 
#%% APERTURE SPACE:  SPHERICAL CONVERGING BEAM
xS = 0; yS = 0; zS = 0.2   # source point S  [m]
a =  0.01 #7.00e-4    #0.01       # radius of aperture  [m]

f = zS                     # focal length                 
zQ = 0
xQ = np.linspace(-a, a, nQ)
yQ = np.linspace(-a, a, nQ)  
XQ,YQ = np.meshgrid(xQ,yQ) 
RQ = (XQ**2 + YQ**2)**0.5

rS = sqrt(xS**2 + yS**2 + zS**2)
rQS = sqrt((XQ-xS)**2 + (YQ-yS)**2 + (zQ-zS)**2) 
EQ = exp(-ik*rQS)/rQS
EQ[RQ > a] = 0

# Comment / uncomment for different aperture fields
# Half-circular aperture
#EQ[YQ < -0] = 0 

# Annular aperture  aI = inner radius
aI = 0.95*a             
EQ[RQ < aI] = 0


IQ = real(np.real(np.conj(EQ)*EQ))   #   IQ(x,y,0)
IQ = IQ/amax(IQ)
IQx = IQ[NQ,:]                       #   IQ(x,y=0,z=0)
IQy = IQ[:,NQ]                       #   IQ(x = 0, y,z=0)
         
#%% OBSERVATION SPACE   [distances m]
u1 = -200.1; u2 = 200.2
uP = linspace(u1,u2,nP)
du = uP[1] - uP[0]
K = zS**2/(k*a**2)
zP = zS + K*uP
xP = 0; yP = 0     
EP = np.zeros(nP)+np.zeros(nP)*1j   
   
# Numerical aperture / Fresnel number
NA = a/sqrt(a**2 + f**2)
NF = a**2/(wL*zS)

#%% SETUP 

# Simpson's [2D] coefficients
S = np.ones(nQ)
R = np.arange(1,nQ,2);   S[R] = 4;
R = np.arange(2,nQ-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy

#%% COMPUATION OF DIFFRACTION INTEGRAL FOR ELECTRIC FIELD and IRRADIANCE
for q in range(nP):
       rPQ = np.sqrt((xP - XQ)**2 + (yP - YQ)**2 + zP[q]**2)
       rPQ3 = rPQ*rPQ*rPQ
       kk = ik * rPQ
       MP1 = exp(kk)
       MP1 = MP1 / rPQ3
       MP2 = zP[q] * (ik * rPQ - 1)
       MP = MP1 * MP2
       EP[q] = sum(sum(EQ*MP*S))

# Intensity (irradiance) optical axis [a.u.]
IPz = real(np.real(EP*np.conj(EP)))
Imax = amax(IPz)
IPz = IPz/Imax
IPzdB = 10*np.log10(IPz)

# Find uP for maximum rradiance 
peak_uP = uP[IPzdB==0][0]

# Find uP for minimums in irradiance 
q = find_peaks(-IPzdB)
dark = uP[q[0]]
dark = dark[dark>0]

#%%  GRAPHICS
# 1   Aperture IQ(x,y=0)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(xQ/a, IQx,'b',lw = 3, label = 'x$_Q$')
ax.plot(yQ/a, IQy,'r',lw = 1, label = 'y$_Q$')
#ax.set_title('m = %0.0f     '%m + 'n = %0.0f' %n,fontsize = 10)
ax.set_xlabel('$x_Q$/a   $y_Q$/a', fontsize=12)
ax.set_ylabel('I$_{Q}$  [a.u.]', fontsize=12)
ax.set_ylim([0,1.01])
#ax.set_xticks(arange(-1,1.2,0.5))
ax.grid()
ax.legend()
fig1.tight_layout()

#%% 2  [2D]  Observation space  IXY(XP, YP, zP)
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3.6,3.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)

cf = ax.pcolor(XQ/a,YQ/a,IQ**sf, cmap='hot')       #Reds
#cf = ax.pcolor(XQ/a,YQ/a,IQ**sf, cmap='Reds')
#cf = ax.contour(VX,VY,IXY**sf, 20,cmap='hot')
fig2.colorbar(cf, ax=ax, location = 'bottom') 
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('$x_{Q}$ / a ', fontsize=10)
ax.set_ylabel('$y_{Q}$  /a ', fontsize=10)
ax.set_xticks(arange(-1,1.2,0.5))
ax.set_yticks(arange(-1,1.2,0.5))

fig2.tight_layout()

#%% 3  GRAPHICS  [1D] Observation space IPz (xP = 0, yP = 0)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(uP, IPz,'b',lw = 3,label = 'RS1')
#ax.plot(u, ID,'r',lw = 1,label = 'Debye')

ax.set_xlabel('$u_P$ ', fontsize=12)
ax.set_ylabel('I$_{Pz}$  [a.u.]', fontsize=12)
ax.grid()
#ax.legend(fontsize = 10)
fig3.tight_layout()

#%% 4   [1D] Observation space IPx(xP,y=0,zP)   dB
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig4, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(uP, IPzdB,'b',lw = 2,label = 'RS1')
#ax.plot(u, IDdB,'r',lw = 1, label = 'Debye')

ax.set_xlabel('$u_P$ ', fontsize=12)
ax.set_ylabel('I$_{Pz}$  [dB]', fontsize=12)
ax.grid()
#ax.legend(fontsize = 10)
fig4.tight_layout()

#%% Console summary
print('  ')
print('emRSFBZ3.py')
print('NQ = %0.0f   '%NQ + 'NP = %0.0f' %NP)
print('nQ = %0.0f   '%nQ + 'nP = %0.0f' %nP)
q = wL*1e9; print('wavelength  wL = %0.0f  nm' %q )
q = a*1e3; print('aperture outside radius  a  = %0.3f  mm' % q ) 
q = aI*1e3; print('aperture inside radius   aI = %0.3f  mm' % q )
print('Source')
print('  xS = %0.3f m  '%xS + 'yS = %0.3f m  ' %yS + 'zS = %0.3f m' %zS)
print('  Focal length f = %0.3f m' %f)
print('Numerical aperture NA = %0.4f   '% NA ) 
print('Fresnel number NF = %0.3f   '% NF ) 
print('Imax = %0.2e  a.u.  '% Imax) 
print('u1 =  %0.2f ' % u1 + 'u2 =  %0.2f' % u2)
print('Uncertainty in uP, du = %0.2f' % du)
print('Max irradiance at uP = %0.4f' %peak_uP)
print('Min at uP')
for q in range(len(dark)):
     print('   %0.2f' % dark[q], end=' ')
tExe = time.time() - tStart
print('\nExecution time %0.0f s' %tExe)

#%%

"""
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')

"""

