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

#%%
import numpy as np
from numpy import pi, exp, linspace, zeros, amax, sqrt, log10 
from numpy.linalg import eig
from scipy.integrate import odeint, quad, dblquad, simps
#from scipy import sqrt
from scipy.special import j1
from scipy.special import jv
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
NQ = 120
NP = 120

nQ = 2*NQ + 1; nP = 2*NP+1 
         
# Wavelength [m]
wL = 632.8e-9 

cL = 2.99792458e8     # speed of light
eps0 = 8.854187e-12   # % permittivity of free space


#%% Aperautre space
a = 4e-4
zQ = 0
xQ = linspace(-a,a,nQ)    
yQ = linspace(-a,a,nQ)
XQ, YQ = np.meshgrid(xQ,yQ)
RQ = (XQ**2 + YQ**2)**0.5
I0 = 1
E0 = sqrt(2*I0/(cL*eps0))
EQ = E0*np.ones([nQ,nQ])
EQ[RQ > a] = 0
IQ = (cL*eps0/2)*np.real(EQ*np.conj(EQ))
IQx = IQ[NQ,:]


#%% Observation spacehalf: half-width % zP   [m] 
zP = 1.1
L = 8*a   
xP = linspace(-L,L,nP)     
yP = linspace(-L,L,nP)
XP, YP = np.meshgrid(xP,yP)
EP = np.zeros([nP,nP])+np.zeros([nP,nP])*1j    


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
         

#%% Irradiance) XY plane, X and Y axes [a.u.] / power
IXY = (cL*eps0/2)*np.real(EP*np.conj(EP))
Iy = IXY[:,NP]
Ix = IXY[NP,:]

IdB = 10*log10(Ix/max(Ix))

#  Power enclosed with a circle [a.u.]
r  = xP[NP:nP]
Ir  = Ix[NP:nP]
Pr = zeros(len(r))

for c in range(len(r)):
    if c > 1:
        Pr[c] = simps(r[0:c]*Ir[0:c],r[0:c])

Pr = 100*Pr/max(Pr)

#%%    Theoretical irradiance IT
sinT = xP/(xP**2 + zP**2)**0.5
q = (2*pi/wL)*sinT
v = q*a+1e-16
IT = (jv(1,v)/(v))**2
ITdB = 10*log10(IT/max(IT))


#%%  Find xP for first dark ring in IXY and xP for minima
q = find_peaks(-IdB)
xmin = xP[q[0]]/a
xDark1 = xmin[xmin>0][0]
PDark1 = Pr[r/a>xDark1][0]
xDark = xmin[xmin>0]

# Relative intensity of maxima
q = find_peaks(Ix)
peaks = Ix[q[0]]/max(Ix)

# Bessel function of 1st kind   rho == v
# Radial optical coordinates for zeros in Bessel function
num = 9999
v = linspace(0,25,num)
J1 = j1(v)

J1index = zeros(10); p = 0
for c in range(num-2):
    q = J1[c]*J1[c+1]
    if q <= 0:
        J1index[p] = c
        p = int(p+1)     
J1index = J1index.astype(int)
vZeros = v[J1index]    

# Rayleigh length  [m]
RL = 4*a**2/wL  
              

#%%  GRAPHICS
# 1  Aperture [1D] radial irradiance
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('$r_Q$ / a ', fontsize=12)
ax.set_ylabel('I$_Q$  ', fontsize=12)
ax.grid()
ax.set_xlim([-1,1]); ax.set_ylim([0,1.2])
ax.plot(xQ/a,IQx,'k',lw = 2)
fig1.tight_layout()

# 2 Aperture [2D] irradiance 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
plt.pcolor(XQ/a,YQ/a,IQ, cmap='Reds')      # Greys_r
ax.set_aspect('equal', adjustable='box')
ax.set_xticks([]); ax.set_yticks([])
fig2.tight_layout()

# 3 Aperture [3D] irradiance
plt.rcParams['font.size'] = 12 
plt.rcParams["figure.figsize"] = (3,3)
fig3 = plt.figure()
ax = fig3.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([]); ax.set_xticks([]); ax.set_yticks([])
ax.set(zlabel=None)
ax.plot_surface(XQ/a,YQ/a,IQ,cmap = 'Reds')
ax.view_init(elev=42, azim=-44, roll=0)
fig3.tight_layout()

#%% 4  Observation space [1D] irradiance Ix
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('$r_P$ / a ', fontsize=12)
ax.set_ylabel('I$_P$  ', fontsize=12)
ax.set_title('z$_p$ = %0.2f  m' %zP)
ax.grid()
ax.plot(xP/a,Ix/max(Ix),'k',lw = 2)
ax.plot(xP/a,IT/max(IT),'r',lw = 1)
fig4.tight_layout()

#%% 5   Observation space [1D] irradiance IxdB
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('$r_P$ / a ', fontsize=12)
ax.set_ylabel('I$_P$  ', fontsize=12)
ax.set_title('z$_p$ = %0.2f  m' %zP)
ax.grid()
ax.plot(xP/a,IdB,'k',lw = 2)
ax.plot(xP/a,ITdB,'r',lw = 1)
fig5.tight_layout()

#%%   6  Observation space [2D] irradiance
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig6, ax = plt.subplots(nrows=1, ncols=1)
plt.pcolor(XP/a,YP/a,IXY**0.23, cmap='Reds')
ax.set_xlabel('$x_P$/a', fontsize=12)
ax.set_ylabel('$y_P$/a', fontsize=12)
ax.set_title('z$_p$ = %0.2f  m' %zP)
ax.set_aspect('equal', adjustable='box')
fig6.tight_layout()

#%%  7  Observation space [3D] irradiance
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3,3)
fig7 = plt.figure()
fig7.subplots_adjust(top = 0.95, bottom = 0.143, left = 0.03,\
                     right = 0.95, hspace = 0.20,wspace = 0.20)
ax = fig7.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([])
ax.set_xlabel('$x_P$/a', fontsize=10)
ax.set_ylabel('$y_P$/a', fontsize=10)
ax.set_title('z$_p$ = %0.2f  m' %zP)
ax.set(zlabel=None)
ax.plot_surface(XP/a,YP/a,IXY**0.2,cmap = 'jet')
ax.view_init(elev=50, azim=-64, roll=0)

#%%   8 power enclosed within rings
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig8, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(r/a,Pr,'k',lw = 2)
ax.set_xlabel('$r_P$ / a ', fontsize=12)
ax.set_ylabel('percentage P  ', fontsize=12)
ax.set_title('z$_p$ = %0.2f  m' %zP)
ax.grid()
ax.plot([xDark1,xDark1],[0,PDark1],'r',lw = 2)
ax.plot([0,xDark1],[PDark1,PDark1],'r',lw = 2)
fig8.tight_layout()

#%%   9  Bessel function of the first kind
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig9, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(v,J1,'k',lw = 2)
ax.set_xlabel('$v_P$ ', fontsize=12)
ax.set_ylabel('J$_1$  ', fontsize=12)
ax.grid()
fig9.tight_layout()
    

#%% Console summary
print('  ')
print('nQ = %0.0f   ' %nQ  + 'nP = %0.0f'%nP)
q = a*1e3; print('aperture radius = %0.3f  mm  ' %q)
q = wL*1e9; print('Wavelength wL = %0.1f  nm  ' %q)
q = L/a; print('Observation space: max(rP/a) = %0.0f ' % q)
print('Observation space: zP = %0.3f m ' % zP)
print('Rayleigh length RL = %0.3f  m' %RL)
print('First dark ring XDark1/a = %0.3f' % xDark1)
print('First dark ring percentage power enclosed  PDark1/a = %0.0f' % PDark1)
print('')
print('Fraunhoffer diffraction: zero irradiance at rP/a')
for c in range(len(xDark)):
     print('   %0.2f' % xDark[c], end = '  ')
print(' '); print('  ')
print('Fraunhoffer diffraction: relative irradiance of the maxima')
for c in range(len(peaks)):
    print(' %0.4f' %peaks[c], end = ' ')
print('  '); print('  ')
  

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

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
fig9.savefig('a9.png')

"""