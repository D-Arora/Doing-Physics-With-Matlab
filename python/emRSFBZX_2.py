# -*- coding: utf-8 -*-
"""
emRSFBZX.py               11May 2025
 
 FOCUSED BEAM  focal length zS
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   IZX irradiance
   x,y,z   [1D] variables   /   X,Y [2D] variables 
   
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRSFB01.pdf
"""

import numpy as np
from numpy import pi, arange, exp, linspace, zeros, amax, sqrt, real
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
NQ = 20
NP = 20

nQ = 2*NQ + 1; nP = 2*NP+1 

sf = 0.25          # scaling factor for IXY plots

# Wavelength [m]
wL = 500e-9  

k = 2*pi/wL            # propagation constant
ik = k*1j              # jk
 
#%% APERTURE SPACE:  SPHERICAL CONVERGING BEAM
xS = 0; yS = 0; zS = 0.2   # source point S  [m]
a = 7.00e-4    #0.01       # radius of aperture  [m]
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
EQ[YQ < -0] = 0
IQ = real(np.real(np.conj(EQ)*EQ))   #   IQ(x,y,0)
IQ = IQ/amax(IQ)
IQx = IQ[NQ,:]                       #   IQ(x,y=0,z=0)
IQy = IQ[:,NQ]                       #   IQ(x = 0, y,z=0)

#%% OBSERVATION SPACE   [distances m]

#u1 = -25.2; u2 = 30.2
#v1 = 0.002; v2 = 30.2
#v1 = -25.2; v2 = -0.002
#v1 = -25.2; v2 = 25.2

yP = 0; YP = 0   
x2 = 1e-4; x1 = -x2
z1 = 0.5*f; z2 = 2*f
xP = linspace(x1,x2,nP)
zP = linspace(z1,z2,nP)
# uP = linspace(u1,u2,nP)
# K = zS**2/(k*a**2)
# zP = zS + K*uP
# vP = linspace(v1,v2,nP)
# K = zS/(k*a)
# xP = K*vP

ZP, XP = np.meshgrid(zP,xP) 
#UP, VP = np.meshgrid(uP,vP) 
EP = np.zeros([nP,nP])+np.zeros([nP,nP])*1j   

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
for c1 in range(nP):
    for c2 in range(nP):
         rPQ = np.sqrt( (XP[c1,c2] - XQ)**2 + (YP - YQ)**2 + ZP[c1,c2]**2 )
         rPQ3 = rPQ*rPQ*rPQ
         kk = ik * rPQ
         MP1 = exp(kk)
         MP1 = MP1 / rPQ3
         MP2 = ZP[c1,c2] * (ik * rPQ - 1)
         MP = MP1 * MP2
         EP[c1,c2] = sum(sum(EQ*MP*S))

# Intensity (irradiance) optical axis [a.u.]
IZX = real(np.real(EP*np.conj(EP)))
Imax = amax(IZX)
IZX = IZX/Imax

Ix = IZX[:,NP]
Iz = IZX[NP,:]

#%%  GRAPHICS
# 1   Aperture IQ(x,y=0)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(xQ/a, IQx,'k',lw = 2)

#ax.set_title('m = %0.0f     '%m + 'n = %0.0f' %n,fontsize = 10)
ax.set_xlabel('$x_Q$/a ', fontsize=12)
ax.set_ylabel('I$_{Qx}$  [a.u.]', fontsize=12)
ax.set_ylim([0.98,1.02])
ax.set_xticks(arange(-1,1,0.5))
ax.grid()
fig1.tight_layout()

#%% 2  [2D]  Observation space  IXY(XP, YP, zP)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)
fig2, ax = plt.subplots(nrows=1, ncols=1)

#cf = ax.contourf(UP,VP,IZX**sf, 12,cmap='hot')
#cf = ax.contour(UP,VP,IZX**sf, 10,cmap='hot')
#cf = ax.pcolormesh(UP,VP,IZX**sf, cmap='hot') 
cf = ax.pcolormesh(ZP,XP*1e6,IZX**sf, cmap='hot')  
fig2.colorbar(cf, ax=ax, location = 'bottom') 

ax.set_xlabel('$z_P$ [m]', fontsize=12)
ax.set_ylabel(r'$r_P$ [$\mu $m]', fontsize=12)
#ax.set_xticks(arange(-10,12,5))
#ax.set_yticks(arange(-10,12,5))
fig2.tight_layout()

#%% 3  [3D] Observation space  IXY(XP, YP, zP)  
# plt.rcParams["figure.figsize"] = (6,4)
# fig3 = plt.figure()
# ax = fig3.add_subplot(1, 1, 1, projection='3d')
# #ax.set_zticks([]); ax.set_xticks([])
# ax.set_zticks([])
# ax.set(zlabel=None)
# #ax.set_title('f = %0.3f     '% f + 'I$_{max}$ = %0.0f' %Imax, fontsize = 12)
# ax.set_xlabel('u$_P$'); ax.set_ylabel('v$_P$')
# ax.plot_surface(UP,VP,IZX**sf,cmap = 'jet')

# ax.view_init(elev=33, azim=-58, roll=0)
# fig3.tight_layout()


#%% Console summary
print('emFBZX.py  ')
print('nQ = %0.0f   '%nQ + 'nP = %0.0f' %nP)
q = wL*1e9; print('wavelength  wL = %0.0f  nm' %q )
q = a*1e3; print('aperture radius  a = %0.3f  mm' % q ) 
print('Source')
print('  xS = %0.2f m  '%xS + 'yS = %0.2f m  ' %yS + 'zS = %0.2f m' %zS)
print('  Focal length f = %0.2f m' %f)
print('Numerical aperture N.A. = %0.5f   '% NA ) 
print('Fresnel number NF = %0.0f   '% NF ) 
#print('max(xP) =  %0.2e m' %Pmax)
#print('zP = %0.3f m   '% zP ) 
print('Imax = %0.2e  a.u.  '% Imax) 
print('Scaling factor  sf = %0.2f ' % sf)

tExe = time.time() - tStart
print('  ')
print('Execution time %0.0f s' %tExe)


#%%

"""

fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')

"""

