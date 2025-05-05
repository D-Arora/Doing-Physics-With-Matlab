# -*- coding: utf-8 -*-
"""
emRSHGXY.py               1 May 2025
 
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
from numpy import pi, arange, exp, linspace, zeros, amax, sqrt 
import matplotlib.pyplot as plt
import time
from scipy import special
from mpl_toolkits.mplot3d import Axes3D
tStart = time.time()

plt.close('all')

#%% INPUT PARAMETERS 
# Grid points: Q aperture space   /   P observation space 
NQ = 60
NP = 60

nQ = 2*NQ + 1; nP = 2*NP+1 

f = 1          # scaling factor for IXY plots

# Wavelength [m]
wL = 632.8e-9  
 
#%% APERTURE SPACE:  HERMITRE-GAUSSS [2D] FUNCTION
m = 2              # x and y orders
n = 2            
a = 4e-4          # radius of circular aperture [m]
w0 = 0.25*a       # standard deviation of Gausssian function
A = 1e-10         # aperture electric field amplitude
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
EQ = A*HH
IQ = np.real(np.conj(EQ)*EQ)         #   IQ(x,y,0)
#IQ = IQ/amax(IQ)
IQx = IQ[NQ,:]                       #   IQ(x,y=0,0)

#%% OBSERVATION SPACE   [distances m]
L = 5

zP = 0.5   # <<<<< 

xPmax = L*a; yPmax = L*a;  
xP = linspace(-xPmax,xPmax,nP)     
yP = linspace(-yPmax,yPmax,nP)
XP, YP = np.meshgrid(xP,yP)    
EP = np.zeros([nP,nP])+np.zeros([nP,nP])*1j         

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
    for c2 in range(nP):
         rPQ = np.sqrt((XP[c1,c2] - XQ)**2 + (YP[c1,c2] - YQ)**2 + zP**2)
         rPQ3 = rPQ*rPQ*rPQ
         kk = ik * rPQ
         MP1 = exp(kk)
         MP1 = MP1 / rPQ3
         MP2 = zP * (ik * rPQ - 1)
         MP = MP1 * MP2
         EP[c1,c2] = sum(sum(EQ*MP*S))

# Intensity (irradiance) optical axis [a.u.]
IXY = np.real(EP*np.conj(EP))
Imax = amax(IXY)
#IXY = IXY/Imax

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

# 2   Aperture  IQ(XQ,YQ)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3.2,3.2)
fig2, ax = plt.subplots(nrows=1, ncols=1)

cf = ax.pcolor(XQ/a,YQ/a,IQ**f, cmap='jet')  

fig2.colorbar(cf, ax=ax, location = 'bottom',shrink=1)    
ax.set_aspect('equal', adjustable='box')
ax.set_title('m = %0.0f     '%m + 'n = %0.0f' %n,fontsize = 10)
ax.set_xticks(arange(-1,1.2,0.5)); ax.set_yticks(arange(-1,1.2,0.5))
fig2.tight_layout()

# 3   Aperture  IQ(XQ,YQ)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3,3)
fig3 = plt.figure()
ax = fig3.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([]); ax.set_xticks([]); ax.set_yticks([])
ax.set(zlabel=None)

ax.plot_surface(XQ/a,YQ/a,IQ**f,cmap = 'jet')

ax.view_init(elev=42, azim=-44, roll=0)
fig3.tight_layout()


#%% 4   Observation space  IXY(XP, YP, zP)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3.2,3.2)
fig4, ax = plt.subplots(nrows=1, ncols=1)

cf = ax.pcolor(XP/a,YP/a,IXY**f, cmap='jet')  

fig4.colorbar(cf, ax=ax, location = 'bottom') 
ax.set_aspect('equal', adjustable='box')
ax.set_title('m = %0.0f     '%m + 'n = %0.0f' %n +
             '     z$_P$ = %0.3f  m' %zP,fontsize = 12)
ax.set_xticks(arange(-5,6,2.5))
ax.set_yticks(arange(-5,6,2.5))
fig4.tight_layout()

# 5 Observation space  IXY(XP, YP, zP)  
plt.rcParams["figure.figsize"] = (3.2,3.2)
fig5 = plt.figure()
ax = fig5.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([]); ax.set_xticks([]); ax.set_yticks([])
ax.set(zlabel=None)
ax.set_title('f = %0.3f     '% f + 'I$_{max}$ = %0.0f' %Imax, fontsize = 12)

ax.plot_surface(XP/a,YP/a,IXY**f,cmap = 'jet')

ax.view_init(elev=42, azim=-44, roll=0)
fig5.tight_layout()


#%% Console summary
print('  ')
print('nQ = %0.0f   '%nQ + 'nP = %0.0f' %nP)
q = wL*1e9; print('wavelength  wL = %0.0f  nm' %q )
q = a*1e3; print('aperture radius  a = %0.2f  mm' % q ) 
print('scaling factor for IXY plots  f = %0.2f   '% f ) 
print('m = %0.0f   '%m + 'n = %0.0f' %n) 
print('zP = %0.3f m   '% zP ) 
print('Imax = %0.2f   '% Imax) 
tExe = time.time() - tStart
print('Execution time %0.0f s' %tExe)


#%%

fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')

"""

"""

