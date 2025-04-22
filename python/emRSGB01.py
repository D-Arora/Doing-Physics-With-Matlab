# -*- coding: utf-8 -*-
"""
emRS02.py               April 2025
 

COMPUTATIONAL OPTICS
   GAUSSIAN BEAM
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   XY irradiance
   x,y,z   [1D] variables   /   X,Y [2D] variables 
   

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
 Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRSGB01.pdf
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
# Gaussian beam s^2 variance
s = 0.5*a
 
# Observation spacehalf: half-width % zP   [m] 
xPmax = 4*a; yPmax = 4*a; zP = 0.05   

#%% SETUP 
k = 2*pi/wL            # propagation constant
ik = k*1j              # jk

# Initialise matrices
unit = np.ones([nQ,nQ])    # unit matrix
rPQ = np.zeros([nQ,nQ]); rPQ3 = np.zeros([nQ,nQ])
MP1 = np.zeros([nQ,nQ]); MP2 = np.zeros([nQ,nQ]); kk = np.zeros([nQ,nQ]) 
MP = np.zeros([nQ,nQ])
EQ = np.ones([nQ,nQ])

#%% Aperture space: Gaussian beam
xQmin = -aQx/2;  xQmax = aQx/2
yQmin = -aQy/2;  yQmax = aQy/2
xQ = linspace(xQmin,xQmax,nQ)    
yQ = linspace(yQmin,yQmax,nQ)
XQ, YQ = np.meshgrid(xQ,yQ)

RQ = (XQ**2 + YQ**2)**0.5

# Aperture electric field and intensity  [a.u.]
EQ = exp(-RQ**2/(2*s**2))
EQ[RQ > xQmax] = 0
IQ = np.real(np.conj(EQ)*EQ)

# Aperture waist  wQ = w0/a    calculated from IQ
xZero = int(np.floor(nQ/2))
IQx = IQ[:,xZero]
wQ = abs(xQ[IQx>exp(-2)][0]/a)

#%% Observation space
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
#IXY = IXY/amax(amax(IXY))
indexXY = num*2   
Iy = IXY[:,indexXY]
Ix = IXY[indexXY,:]
#Ix = Ix/max(Ix); Iy = Iy/max(Iy)
Iref = 1.47337143616115e24    # only for nQ = 199  nP = 59
Ix = Ix/Iref; Iy = Iy/Iref
# Observation space waist
wP = abs(xP[Ix>max(Ix)*exp(-2)][0]/a)

# Waist and beam spot:theory
w0 = wQ*a
zR = pi*w0**2/wL

wT = w0*(1+(zP/zR)**2)**0.5/a


#%%  GRAPHICS
# 1  Aperture intensity 
plt.rcParams["figure.figsize"] = (5,4)
fig1 = plt.figure()
ax = fig1.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([]); ax.set_xticks([]); ax.set_yticks([])
ax.set(zlabel=None)
ax.set_title('Gaussian beam')
ax.plot_surface(XQ/a,YQ/a,IQ,cmap = 'jet')
ax.view_init(elev=42, azim=-44, roll=0)
fig1.tight_layout()

#%%  2 Aperture intensity 
plt.rcParams["figure.figsize"] = (4,4)
fig2, ax = plt.subplots(nrows=1, ncols=1)
plt.pcolor(XQ/a,YQ/a,IQ, cmap='Greys_r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('$x_Q$ / a ', fontsize=12)
ax.set_ylabel('$y_Q$ / a ', fontsize=12)
fig2.tight_layout()

#%%   3 aperture radial intensity
xZero = int(np.floor(nQ/2))
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(xQ/a,IQ[xZero,:],'k',lw = 2)
ax.plot([0,wQ],[exp(-2),exp(-2)],'r',lw=1)
ax.plot([wQ,wQ],[0,exp(-2)],'r',lw=1)
q = wQ; ax.set_title('waist   w$_0$/a = %0.3f ' % q,fontsize = 12)
ax.set_xlabel('$r_Q$ / a ', fontsize=12)
ax.set_ylabel('I$_Q$  [a.u.]', fontsize=12)
ax.grid()
fig3.tight_layout()

#%%
# 4  Observation space intensity 
plt.rcParams["figure.figsize"] = (5,4)
fig4 = plt.figure()
ax = fig4.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([]); ax.set_xticks([]); ax.set_yticks([])
ax.set(zlabel=None)
ax.plot_surface(XP/a,YP/a,IXY,cmap = 'Greys')
ax.view_init(elev=42, azim=-44, roll=0)
fig4.tight_layout()

#%%  5 Observation space intensity 
plt.rcParams["figure.figsize"] = (2,2)
fig5, ax = plt.subplots(nrows=1, ncols=1)
plt.pcolor(XP/a,YP/a,IXY, cmap='Greys_r')
ax.set_aspect('equal', adjustable='box')
#ax.set_xlabel('$x_P$ / a ', fontsize=12)
#ax.set_ylabel('$y_P$ / a ', fontsize=12)
#ax.set_xticks(np.arange(-4,5,2))
#ax.set_yticks(np.arange(-4,5,2)) 
ax.set_xticks([]);ax.set_yticks([]) 
ax.text(-3.6,-3.6,'zP = %0.2f m' %zP,color = [1,1,0],fontsize = 10)            
fig5.tight_layout()

#%%   6 observation space intensity
xZero = int(np.floor(nQ/2))
plt.rcParams["figure.figsize"] = (6,3)
fig6, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(xP/a,Ix,'k',lw = 2)
ax.plot([0,wP],[max(Ix)*exp(-2),max(Ix)*exp(-2)],'r',lw=1)
ax.plot([wP,wP],[0,max(Ix)*exp(-2)],'r',lw=1)

ax.set_title('z$_P$ = %0.2f m    ' %zP
             + 'w$_P$/a = %0.3f   ' % wP 
             + 'w$_T$/a =%0.3f ' %wT,fontsize = 12)

ax.set_xlabel('$r_P$ / a ', fontsize=12)
ax.set_ylabel('I$_P$  [a.u.]', fontsize=12)
ax.grid()
fig6.tight_layout()


#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png') 
fig5.savefig('a5.png')
fig6.savefig('a6.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')


#%% Console summary
print('  ')
print('nQ = %0.0f   ' %nQ  + 'nP = %0.0f'%nP)
q = a*1e3; print('aperture radius = %0.3f  mm  ' %q)
q = s*1e3; print('Gaussian beam: width s = %0.3f  mm  ' %q)
q = wL*1e9; print('Wavelength wL = %0.1f  nm  ' %q)
print('Observation space: zP = %0.3f  mm  ' % zP)
print('  ')
q = w0/a; print('Gaussian beam waist w0/a = %0.3f' %q)
q = wT; print('Gaussian beam spot (theoretical) wT/a = %0.3f' %q)
q = max(Ix);print('Relative max irradiance (normalized to 1 at zP = 1.0 m) = %0.3f' % q)        
print('  ')
print('Execution time %0.0f  s ' %tExe)
