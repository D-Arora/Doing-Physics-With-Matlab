# -*- coding: utf-8 -*-
"""
emRSFBGXY1.py               12 May 2025
 
 FOCUSED BEAM  focal length rS
   BEAM PROPAGATION FROM A PLANAR APERTURE
   RAYLEIGH-SOMMERFELD 1 DIFFRACTION INTEGRAL SOLVED
   CIRCULAR APERTURE  radius a
   IXY irradiance
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
from scipy.integrate import simps
tStart = time.time()

plt.close('all')

#%% INPUT PARAMETERS 
# Grid points: Q aperture space   /   P observation space 
NQ = 100
NP = 100

nQ = 2*NQ + 1; nP = 2*NP+1 

sf = 0.25          # scaling factor for IXY plots

# Wavelength [m]
wL = 500e-9  

k = 2*pi/wL            # propagation constant
ik = k*1j              # jk
 
#%% APERTURE SPACE:  SPHERICAL CONVERGING BEAM
a = 0.01    # radius of aperture  [m]
zQ = 0

xQ = np.linspace(-a, a, nQ)
yQ = xQ
XQ,YQ = np.meshgrid(xQ,yQ) 
RQ = (XQ**2 + YQ**2)**0.5

# SOURCE POINT   [m]
xS = 0; yS = 0; zS = 0.02
# Focal length
f = zS

rS = sqrt(xS**2 + yS**2 + zS**2)
rQS = sqrt((XQ-xS)**2 + (YQ-yS)**2 + (zQ-zS)**2) 
EQ = exp(-ik*rQS)/rQS
EQ[RQ > a] = 0
IQ = real(np.real(np.conj(EQ)*EQ))         #   IQ(x,y,0)
IQ = IQ/amax(IQ)
IQx = IQ[NQ,:]                       #   IQ(x,y=0,0)

#%% OBSERVATION SPACE   [distances m]
zP = zS
v1 = 25
vP = linspace(-v1,v1,nP)
xP = vP*zS/(k*a)
yP = xP
VX, VY = np.meshgrid(vP,vP) 
XP, YP = np.meshgrid(xP,yP)    
EP = np.zeros([nP,nP])+np.zeros([nP,nP])*1j   

# Numerical aperutre andFresnel number
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
         rPQ = np.sqrt((XP[c1,c2] - XQ)**2 + (YP[c1,c2] - YQ)**2 + zP**2)
         rPQ3 = rPQ*rPQ*rPQ
         kk = ik * rPQ
         MP1 = exp(kk)
         MP1 = MP1 / rPQ3
         MP2 = zP * (ik * rPQ - 1)
         MP = MP1 * MP2
         EP[c1,c2] = sum(sum(EQ*MP*S))

# Intensity (irradiance) optical axis [a.u.]
IXY = real(np.real(EP*np.conj(EP)))
Imax = amax(IXY)
IXY = IXY/Imax
IPx = IXY[NP,:]                       #   IP(x,y=0,zS)
IPxdB = 10*np.log10(IPx)

# Bessel function
v = linspace(0.01,20,999)
J1 = jv(1, v)
IB = (J1/v)**2
IB = IB/max(IB)
IBdB = 10*np.log10(IB)

q = find_peaks(-IBdB)
J1zeros = v[q[0]]

#  Power enclosed with a circle [a.u.]
r  = xP[NP:nP]
vr = vP[NP:nP]
Ir  = IPx[NP:nP]
Pr = zeros(len(r))

for c in range(len(r)):
    if c > 1:
        Pr[c] = simps(vr[0:c]*Ir[0:c],vr[0:c])

Pr = 100*Pr/max(Pr)

# Find vP for dark rings in irradiance 
q = find_peaks(-IPxdB)
dark = vP[q[0]]
dark = dark[dark>0]
#dark0 = dark[dark>0][0]
dark0 = 3.77  #  vP(dark) Fraunhoffer first dark ring
Pdark = Pr[vr>dark0][0]


#%%  GRAPHICS
# 1   Aperture IQ(x,y=0)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(xQ/a, IQx,'k',lw = 2)

#ax.set_title('m = %0.0f     '%m + 'n = %0.0f' %n,fontsize = 10)
ax.set_xlabel('$x_Q$/a ', fontsize=12)
ax.set_ylabel('I$_{Qx}$  [a.u.]', fontsize=12)
ax.set_ylim([0.6,1.02])
ax.set_xticks(arange(-1,1.4,0.5))
ax.grid()
fig1.tight_layout()

#%% 2  [2D]  Observation space  IXY(XP, YP, zP)
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3.3,3.3)
fig2, ax = plt.subplots(nrows=1, ncols=1)

cf = ax.pcolor(VX,VY,IXY**sf, cmap='Reds')  

fig2.colorbar(cf, ax=ax, location = 'bottom') 
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('$v_{Px}$ ', fontsize=10)
ax.set_ylabel('$v_{Py}$ ', fontsize=10)
#ax.set_xticks(arange(-10,12,5))
#ax.set_yticks(arange(-10,12,5))
#ax.set_xticks(arange(-20,22,10))
#ax.set_yticks(arange(-20,22,10))

fig2.tight_layout()

#%% 3  [3D] Observation space  IXY(XP, YP, zP)  
plt.rcParams["figure.figsize"] = (4,4)
fig3 = plt.figure()
ax = fig3.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([]); ax.set_xticks([]); ax.set_yticks([])
ax.set(zlabel=None)
#ax.set_title('f = %0.3f     '% f + 'I$_{max}$ = %0.0f' %Imax, fontsize = 12)

ax.plot_surface(VX,VY,IXY**sf,cmap = 'jet')

ax.view_init(elev=50, azim=-43, roll=0)
fig3.tight_layout()

#%% 4   [1D] Observation space IPx(xP,y=0,zP)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig4, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(vP, IPx,'k',lw = 2)


#ax.set_title('m = %0.0f     '%m + 'n = %0.0f' %n,fontsize = 10)
ax.set_xlabel('$v_P$ ', fontsize=12)
ax.set_ylabel('I$_{Px}$  [a.u.]', fontsize=12)
ax.grid()
fig4.tight_layout()

#%% 5   [1D] Observation space IPx(xP,y=0,zP)   dB
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig5, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(vP, IPxdB,'k',lw = 2)
ax.plot(v, IBdB,'b',lw = 2)

ax.set_xlabel('$v_P$ ', fontsize=12)
ax.set_ylabel('I$_{Px}$  [dB]', fontsize=12)
ax.grid()
fig5.tight_layout()

#%%   6 power enclosed within rings
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig6, ax = plt.subplots(nrows=1, ncols=1)

ax.plot(vr,Pr,'k',lw = 2)
ax.plot([dark0,dark0],[0,Pdark],'b',lw = 2)
ax.plot([0,dark0],[Pdark,Pdark],'b',lw = 2)
ax.set_xlabel('$v_P$ ', fontsize=12)
ax.set_ylabel('percentage P  ', fontsize=12)
ax.set_title('vP(dark) = %0.2f' %dark0 +
             '   percent Pdark = %0.1f' %Pdark, color = 'b', fontsize = 10)
ax.grid()

fig6.tight_layout()

#%% Console summary
print('emRSFBXY.py  ')
print('nQ = %0.0f   '%nQ + 'nP = %0.0f' %nP)
q = wL*1e9; print('wavelength  wL = %0.0f  nm' %q )
q = a*1e3; print('aperture radius  a = %0.2f  mm' % q ) 
print('Source')
print('  xS = %0.2f m  '%xS + 'yS = %0.2f m  ' %yS + 'zS = %0.2f m' %zS)
print('  Focal length f = %0.2f m' %f)
print('Numerical aperture N.A. = %0.3f   '% NA ) 
print('Fresnel number NF = %0.0f   '% NF ) 
#print('max(xP) =  %0.2e m' %Pmax)
print('zP = %0.5f m   '% zP ) 
print('Imax = %0.2e  a.u.  '% Imax) 
print('Dark rings vP')
for q in range(len(dark)):
    print('   %0.2f' % dark[q], end=' ')
print('\nDark rings Bessel function zeros')
for q in range(len(J1zeros)):
    print('   %0.2f' % J1zeros[q], end=' ')
print('\nFraunhoffer first dark ring vP(dark) = %0.2f' %dark0) 
#print('Percentage power within first dark ring Pdark = %0.1f' %Pdark) 
   
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
fig6.savefig('a6.png')


"""

