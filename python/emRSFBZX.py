# -*- coding: utf-8 -*-
"""
emRSFBGZX.py               22 May 2025
 
 FOCUSED BEAM  focal length f = zS
   BEAM PROPAGATION FROM A PLANAR CIRCULAR APERTURE radius a
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
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
from scipy.integrate import simps
from scipy.signal import find_peaks
from scipy.special import jv
from mpl_toolkits.mplot3d import Axes3D
tStart = time.time()

plt.close('all')

#%% INPUT PARAMETERS 
# Grid points: Q aperture space   /   P observation space 
NQ = 99; NP = 99
# scaling factor for IXY plots
sf = 0.25          
# radius of aperture  [m]       0.01  0.02   7.00e-4
a =  0.01 #7.00e-4    #0.01 
# source point S  [m]
xS = 0; yS = 0; zS = 1 
# Observation space: Z, X axis limits - optical coordinates
u1 = -60.2; u2 = 60.2
v1 = -40.2; v2 = 40.2

# Wavelength [m]
wL = 500e-9  
# Annular aperture  aI = inner radius
aI = 0.95*a   

#%% SELECT APERTURE FUNCTION FOR EQ using variable A
#   A = 1  --> default value: circular aperture with no mask function applied
#   A = 2  --> half-circular aperture (semi-circular aperture)
#   A = 3  --> annular aperture" inner radius aI
#   A = 4  --> spherical aberration
#   A = 5  --> comatic aberration
A = 1

#   AP = 0  --> do not plot Fraunhoffer diffraction patterns
#   AP = 1  --> add Fraunhoffer diffraction patterns to plots
AP = 0 

# [2D] aperture plot
# QP = 0 --> plot aperture irradiance IQ
# QP = 1 --> plot phase of apertute electric field EQ
QP = 1 

#%% SETUP
nQ = 2*NQ + 1; nP = 2*NP+1     # Grid
k = 2*pi/wL; ik = k*1j         # Propagation constant / jk
f = zS                         # Focal length     
# Numerical aperture / Fresnel number
NA = a/sqrt(a**2 + f**2); NF = a**2/(wL*zS)

#%% APERTURE SPACE:  SPHERICAL CONVERGING BEAM  
zQ = 0
xQ = np.linspace(-a-1e-9, a, nQ)
yQ = np.linspace(-a-1e-9, a, nQ)  
XQ,YQ = np.meshgrid(xQ,yQ) 
RQ = (XQ**2 + YQ**2)**0.5

rS = sqrt(xS**2 + yS**2 + zS**2)
rQS = sqrt((XQ-xS)**2 + (YQ-yS)**2 + (zQ-zS)**2) 
EQ = exp(-ik*rQS)/rQS

# APERTURE MASK FUNCTION

if A == 2:    # Half-circular aperture
   EQ[YQ < -0] = 0 

if A == 3:   # Annular aperture  aI = inner radius
   EQ[RQ < aI] = 0

if A == 4:   # Spherical aberration
   phi = -pi*RQ**4
   T = exp(ik*phi)
   EQ = T*EQ

if A == 5: # COMA aberration
   cosT = YQ/(RQ)
   phi = 6*pi*RQ**3*cosT
   T = exp(ik*phi)
   EQ = T*EQ

EQ[RQ > a] = 0
IQ = real(np.real(np.conj(EQ)*EQ))   #   IQ(x,y,0)
IQ = IQ/amax(IQ)
IQx = IQ[NQ,:]                       #   IQ(x,y=0,z=0)
IQy = IQ[:,NQ]                       #   IQ(x=0, y,z=0)
IQdB = 10*np.log10(IQ+1e-16)         #   decibels
phase = np.angle(EQ)                 #   electric field phase
         


#%% OBSERVATION SPACE   [distances m]
yP = 0; YP = 0   
uP = linspace(u1,u2,nP)
K = zS**2/(k*a**2)
zP = zS + K*uP
vP = linspace(v1,v2,nP)
K = zS/(k*a)
xP = K*vP

ZP, XP = np.meshgrid(zP,xP) 
UP, VP = np.meshgrid(uP,vP) 
EP = np.zeros([nP,nP])+np.zeros([nP,nP])*1j   

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
Ix = IZX[:,NP]; Iz = IZX[NP,:]
IxdB = 10*np.log10(Ix); IzdB = 10*np.log10(Iz)

#%%  GRAPHICS
# 1  [1D] Aperture IQ(x,y,0)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(xQ/a, IQx,'b',lw = 3, label = 'x$_Q$')
ax.plot(yQ/a, IQy,'r',lw = 1, label = 'y$_Q$')
ax.set_xlabel('$x_Q$/a   $y_Q$/a', fontsize=12)
ax.set_ylabel('I$_{Q}$  [a.u.]', fontsize=12)
ax.set_ylim([0,1.01])
#ax.set_xticks(arange(-1,1.2,0.5))
ax.grid()
ax.legend(ncols = 2)
fig1.tight_layout()

#%% 2  [2D]  Aperture space  
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3.6,3.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)


if QP == 0: # Irradiance IQ  [dB]
   cf = ax.pcolor(XQ/a,YQ/a,IQdB, cmap='jet')       
if QP == 1:   # electric field phase [rad/pi]
   cf = ax.pcolor(XQ/a,YQ/a,phase/pi, cmap='jet') 

fig2.colorbar(cf, ax=ax, location = 'bottom') 
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('$x_{Q}$ / a ', fontsize=10)
ax.set_ylabel('$y_{Q}$  /a ', fontsize=10)
ax.set_xticks(arange(-1,1.2,0.5))
ax.set_yticks(arange(-1,1.2,0.5))
fig2.tight_layout()

#%% 3  [2D]  Observation space  IXY(XP, YP, zP)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)

#cf = ax.contourf(UP,VP,IZX**sf, 32,cmap='jet')

#cf = ax.contour(UP,VP,IZX**sf, 20,cmap='jet')
#cf = ax.contour(UP,VP,IZX**sf, 16,colors='k')
cf = ax.pcolormesh(UP,VP,IZX**sf, cmap='jet')  
#cf = ax.pcolormesh(ZP,XP,IZX**sf, cmap='hot')  
fig3.colorbar(cf, ax=ax, location = 'bottom',ticks = arange(0,1.1,0.2)) 

ax.set_xlabel('$u_P$ ', fontsize=12)
ax.set_ylabel('$v_P$ ', fontsize=12)
#ax.set_xticks(arange(-10,12,5))
#ax.set_yticks(arange(-10,12,5))
fig3.tight_layout()

#%% 4  [3D] Observation space  IXY(XP, YP, zP)  
plt.rcParams["figure.figsize"] = (5,3.5)
fig4 = plt.figure()
ax = fig4.add_subplot(1, 1, 1, projection='3d')
#ax.set_zticks([]); ax.set_xticks([])
ax.set_zticks([])
ax.set(zlabel=None)
#ax.set_title('f = %0.3f     '% f + 'I$_{max}$ = %0.0f' %Imax, fontsize = 12)
ax.set_xlabel('u$_P$'); ax.set_ylabel('v$_P$')
ax.plot_surface(UP,VP,IZX**sf,cmap = 'jet')

ax.view_init(elev=33, azim=-58, roll=0)
fig4.tight_layout()

#%% 5  [1D] Observation space IP
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.8,2.5)
fig5, ax = plt.subplots(nrows=1, ncols=2)

ax[0].plot(vP, IxdB,'b',lw = 2)
ax[0].set_xlabel('$v_P$', fontsize=12)
ax[0].set_ylabel('I$_{Px}$  [a.u.]', fontsize=12)
ax[0].grid()

ax[1].plot(uP, IzdB,'r',lw = 2)
ax[1].set_xlabel('$u_P$', fontsize=12)
ax[1].set_ylabel('I$_{Pz}$  [a.u.]', fontsize=12)
ax[1].grid()
fig5.tight_layout()

#%% Console summary
print('emFBZX.py  ')
print('nQ = %0.0f   '%nQ + 'nP = %0.0f' %nP)
q = wL*1e9; print('wavelength  wL = %0.0f  nm' %q )
q = a*1e3; print('aperture radius  a = %0.2f  mm' % q ) 
print('Source')
print('  xS = %0.2f m  '%xS + 'yS = %0.2f m  ' %yS + 'zS = %0.2f m' %zS)
print('  Focal length f = %0.2f m' %f)
print('Numerical aperture N.A. = %0.3f   '% NA ) 
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
fig4.savefig('a4.png')
fig5.savefig('a5.png')


"""

