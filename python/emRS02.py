# -*- coding: utf-8 -*-
"""
emRS02.py               April 2025
 

COMPUTATIONAL OPTICS
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   XY irradiance
   x,y,z   [1D] variables   /   X,Y [2D] variables 
   

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
 Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRS02.pdf
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

tStart = time.time()

#%% INPUT PARAMETERS 
# Grid points: Q aperture space    nQ format odd number
#              P observation space nP format integer * 4 + 1
nQ = 199    
num = 59;           
nP = num*4+1
         
# Wavelength [m]
wL = 632.8e-9 
       
# Aperautre space: radius a / XY dimensions of aperture [m]
a = 4e-4
aQx = 2*a; aQy = 2*a           
 
# Observation spacehalf: half-width % zP   [m] 
xPmax = 10*a; yPmax = 10*a; zP = 1.0       


#%% SETUP 
k = 2*pi/wL            # propagation constant
ik = k*1j              # jk

# Initialise matrices
unit = np.ones([nQ,nQ])    # unit matrix
rPQ = np.zeros([nQ,nQ]); rPQ3 = np.zeros([nQ,nQ])
MP1 = np.zeros([nQ,nQ]); MP2 = np.zeros([nQ,nQ]); kk = np.zeros([nQ,nQ]) 
MP = np.zeros([nQ,nQ])
EQ = np.ones([nQ,nQ])

# Aperture space
xQmin = -aQx/2;  xQmax = aQx/2
yQmin = -aQy/2;  yQmax = aQy/2
xQ = linspace(xQmin,xQmax,nQ)    
yQ = linspace(yQmin,yQmax,nQ)
XQ, YQ = np.meshgrid(xQ,yQ)

RQ = (XQ**2 + YQ**2)**0.5

# Aperture electric field and intensity  [a.u.]
EQ[RQ > xQmax] = 0
IQ = np.real(np.conj(EQ)*EQ)

# Observation space
xPmin = -xPmax; yPmin =  -yPmax
EP = np.zeros([nP,nP])+np.zeros([nP,nP])*1j         

xP = linspace(xPmin,xPmax,nP)     
yP = linspace(yPmin,yPmax,nP)
XP, YP = np.meshgrid(xP,yP)

# Simpson [2D] coefficients
S = np.ones(nQ)
R = np.arange(1,nQ,2);   S[R] = 4;
R = np.arange(2,nQ-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy

#%% COMPUATION OF DIFFRACTION INTEGRAL FOR ELECTRIC FIELD % IRRADIANCE
for c1 in range(nP):
    for c2 in range(nP):
         rPQ = np.sqrt((XP[c1,c2] - XQ)**2 + (YP[c1,c2] - YQ)**2 + zP**2)
         rPQ3 = rPQ*rPQ*rPQ
         kk = ik * rPQ
         MP1 = exp(kk)
         MP1 = MP1 / rPQ3
         MP2 = zP * (ik * rPQ - unit)
         MP = MP1 * MP2
         EP[c1,c2] = sum(sum(EQ*MP*S))

# Intensity (irradiance) XY plane [a.u.]
IXY = np.real(EP*np.conj(EP))
IXY = IXY/amax(amax(IXY))
indexXY = num*2   
Iy = IXY[:,indexXY]
Ix = IXY[indexXY,:]
Ix = Ix/max(Ix); Iy = Iy/max(Iy)

# Power enclosed with a circle [a.u.]
r  = xP[indexXY:nP]
Ir  = Ix[indexXY:nP]
Pr = zeros(len(r))

for c in range(len(r)):
    if c > 1:
       Pr[c] = simps(r[0:c]*Ir[0:c],r[0:c])

Pr = 100*Pr/max(Pr)

#%%
# Observation space: radial positions for zero intensity xZ
#   Zeros for Bessel function of first kind  J1(rho = 0) 
#   Angles (theta) for zeros in diffraction pattern    
rho = np.array([3.8317,7.0156,10.1735,13.3237,16.4706,19.6159])
theta = np.arcsin(rho/(a*k))
xZ = zP*np.tan(theta)/a

# RS1 predictions for location radial positions for zero intensity
#   Radial intenisty array  IX and IdB 
IdB = 10*np.log10(Ix)
q = find_peaks(-IdB)
xZ_RS1 = xP[q[0]]/a
x1 = xZ_RS1[xZ_RS1>0][0]


# Relative intensities of maxima
q = find_peaks(Ix)
peaks = Ix[q[0]]

# Bessel function of 1st kind   rho == v
num = 9999
v = linspace(0,25,num)
J1 = j1(v)

# Radial optical coordinates for zeros in Bessel function
J1index = zeros(10); p = 0
for c in range(num-2):
    q = J1[c]*J1[c+1]
    if q <= 0:
       J1index[p] = c
       p = int(p+1)     
J1index = J1index.astype(int)
vZeros = v[J1index]    

# Rayleigh length  [m]
dRL = 4*a**2/wL 
         
              
#%% Console summary
print('Relative intensities of the maxima in Fraunhoffer diffraction pattern')
for c in range(len(peaks)):
    print('%0.4f' %peaks[c])
print('  ')
print('Zero in intensities in Fraunhoffer diffraction pattern')
for c in range(len(rho)):
    print('%0.2f' % xZ[c])
print(' ')
print('Zero in intensities in Fraunhoffer diffraction pattern')
for c in range(len(xZ_RS1)):
    print('%0.2f' % xZ_RS1[c])
print(' ')
print('Zero in intensities in Fraunhoffer diffraction pattern')
for c in range(7):
    print('%0.3f' % vZeros[c])    
    

#%%  GRAPHICS
# 1  Aperture intensity 
plt.rcParams["figure.figsize"] = (5,4)
fig1 = plt.figure()
ax = fig1.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([]); ax.set_xticks([]); ax.set_yticks([])
ax.set(zlabel=None)
ax.plot_surface(XQ/a,YQ/a,IQ,cmap = 'Greys')
ax.view_init(elev=42, azim=-44, roll=0)
fig1.tight_layout()

#%%  2 Aperture intensity 
plt.rcParams["figure.figsize"] = (3,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
plt.pcolor(XQ/a,YQ/a,IQ, cmap='Greys_r')
ax.set_aspect('equal', adjustable='box')
ax.set_xticks([]); ax.set_yticks([])
fig2.tight_layout()

#%%   3 radial intensiy
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)

fig3, ax = plt.subplots(nrows=1, ncols=2)
ax[0].plot(xP/a,Ix,'k',lw = 2)
ax[0].set_xlabel('$r_P$ / a ', fontsize=12)
ax[0].set_ylabel('I  [a.u.]', fontsize=12)
ax[0].grid()

ax[1].plot(xP/a,10*np.log(Ix),'k',lw = 2)
ax[1].set_xlabel('$r_P$ / a ', fontsize=12)
ax[1].set_ylabel('I  [a.u.]', fontsize=12)
ax[1].grid()
fig3.tight_layout()

#%%   4   pcolor plot  IXY
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
plt.pcolor(XP/a,YP/a,IXY**0.23, cmap='Greys_r')
ax.set_xlabel('$x_P$/a', fontsize=12)
ax.set_ylabel('$y_P$/a', fontsize=12)
ax.set_aspect('equal', adjustable='box')
fig4.tight_layout()

#%%  5  surface plot IXY
plt.rcParams["figure.figsize"] = (5,4)
fig5 = plt.figure()
ax = fig5.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([])
ax.set_xlabel('$x_P$/a', fontsize=10)
ax.set_ylabel('$y_P$/a', fontsize=10)
#ax.set_title('max(S$_{xy}$) = %0.0f  W.m$^{-2}$' %amax(Sxy),fontsize = 10)
ax.set(zlabel=None)
ax.plot_surface(XP/a,YP/a,IXY**0.2,cmap = 'jet')
fig5.tight_layout()

#%%   6 power
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)

fig6, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(r/a,Pr,'k',lw = 2)
ax.set_xlabel('$r_P$ / a ', fontsize=12)
ax.set_ylabel('percentage P  ', fontsize=12)
ax.grid()
ax.set_xlim([0,12]); ax.set_ylim([0,100])

q = Pr[r/a>x1][0]
ax.plot([x1,x1],[0,q],'r',lw = 2)
ax.plot([0,x1],[q,q],'r',lw = 2)

fig6.tight_layout()

#%%   7  Bessel function of the first kind
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)

fig7, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(v,J1,'k',lw = 2)
ax.set_xlabel('$v_P$ ', fontsize=12)
ax.set_ylabel('J$_1$  ', fontsize=12)
ax.grid()
fig7.tight_layout()
    
#%%   8 radial intensiy
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,2.5)
fig8, ax = plt.subplots(nrows=1, ncols=1)
ax.set_title('z$_p$ = %0.2f  m' %zP)
ax.plot(xP/a,Ix,'k',lw = 2)
ax.set_xlabel('$r_P$ / a ', fontsize=12)
ax.set_ylabel('I  [a.u.]', fontsize=12)
ax.grid()
fig8.tight_layout()

#%%
"""
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png') 
fig5.savefig('a5.png')
fig6.savefig('a6.png')
fig7.savefig('a7.png')
fig8.savefig('a8.png')
"""

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
