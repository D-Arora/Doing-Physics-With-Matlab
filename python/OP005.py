# -*- coding: utf-8 -*-
"""
op005.py               Oct 2024
 
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/op005.pdf
"""

import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax 
from numpy.linalg import eig
from scipy.integrate import odeint, quad, dblquad, simps
from scipy import pi, sqrt
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import axes3d
import sympy as sym


tStart = time.time()

# INPUT PARAMETERS -------------------------------------------------------
num = 60             # number for observation space
nP = num*4+1         # observation points for P:  format integer * 4 + 1
nQ = 201              # aperture points for Q  must be ODD

wL = 632.8e-9        # wavelength [m]

# Aperautre space: full-width [m]
aQx = 1e-4; aQy = 2e-4           
 
#Observation spacehalf: half-width % zP   [m] 
xPmax = 2e-2; yPmax = 2e-2; zP = 1      
         

# SETUP -----------------------------------------------------------------
k = 2*pi/wL            # propagation constant
ik = k*1j              # jk

# Initialise matrices
unit = np.ones([nQ,nQ])    # unit matrix
rPQ = np.zeros([nQ,nQ]); rPQ3 = np.zeros([nQ,nQ])
MP1 = np.zeros([nQ,nQ]); MP2 = np.zeros([nQ,nQ]); kk = np.zeros([nQ,nQ]) 
MP = np.zeros([nQ,nQ])

# Aperture space
xQmin = -aQx/2;  xQmax = aQx/2
yQmin = -aQy/2;  yQmax = aQy/2
xQ1 = linspace(xQmin,xQmax,nQ)    
yQ1 = linspace(yQmin,yQmax,nQ)
xQ, yQ = np.meshgrid(xQ1,yQ1)

# Observation space
xPmin = -xPmax; yPmin =  -yPmax
EP = np.zeros([nP,nP])+np.zeros([nP,nP])*1j         

xP1 = linspace(xPmin,xPmax,nP)     
yP1 = linspace(yPmin,yPmax,nP)
xP, yP = np.meshgrid(xP1,yP1)


# Simpson [2D] coefficients
S = np.ones(nQ)
R = np.arange(1,nQ,2);   S[R] = 4;
R = np.arange(2,nQ-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy

#%% COMPUATION OF DIFFRACTION INTEGRAL FOR ELECTRIC FIELD % IRRADIANCE 

for c1 in range(nP):
    for c2 in range(nP):
         rPQ = np.sqrt((xP[c1,c2] - xQ)**2 + (yP[c1,c2] - yQ)**2 + zP**2)
         rPQ3 = rPQ*rPQ*rPQ
         kk = ik * rPQ
         MP1 = exp(kk)
         MP1 = MP1 / rPQ3
         MP2 = zP * (ik * rPQ - unit)
         MP = MP1 * MP2
         EP[c1,c2] = sum(sum(MP*S))

Irr = np.real(EP*np.conj(EP))
Irr = Irr/amax(amax(Irr))
indexXY = num*2+1   
IY = Irr[:,indexXY]
IX = Irr[indexXY,:]

# Position of first zero in irradiance
x0 = wL*zP/(1*aQx)*1e3    # [mm]
y0 = wL*zP/(1*aQy)*1e3    # [m]


#%%  GRAPHICS 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3.9,3)


fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(xP1*1e3,IX,'b',lw = 2, label = 'X')
ax.plot(yP1*1e3,IY,'r',lw = 2, label = 'Y')
ax.set_xlabel('$x_P$ (blue)   $y_P$ (red)  [mm]', fontsize=12)
ax.set_ylabel('I  [a.u.]', fontsize=12)
ax.set_xlim([-20,20])
ax.xaxis.grid()
ax.yaxis.grid()
#ax.legend()
fig.tight_layout()
fig.savefig('a1.png')

fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(xP1*1e3,np.log(IX),'b',lw = 2, label = 'X')
ax.plot(yP1*1e3,np.log(IY),'r',lw = 2, label = 'Y')
ax.set_xlabel('$x_P$ (blue)   $y_P$ (red)  [mm]', fontsize=12)
ax.set_ylabel('I  [dB]', fontsize=12)
ax.set_xlim([-20,20])
ax.xaxis.grid()
ax.yaxis.grid()
#ax.legend()
fig.tight_layout()
fig.savefig('a2.png')


fig, ax = plt.subplots(nrows=1, ncols=1)
#plt.contourf(xP,yP,Irr**0.25, cmap='Reds')
#plt.pcolor(xP,yP,Irr**0.25, cmap='Reds')
plt.pcolor(xP*1e3,yP*1e3,Irr**0.3, cmap='Greys_r')
ax.set_xlabel('$x_P$  [mm]', fontsize=12)
ax.set_ylabel('$y_P$  [mm]', fontsize=12)
ax.set_aspect('equal', adjustable='box')
fig.tight_layout()
fig.savefig('a3.png')


#%%   NUMERICAL SUMMARY
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3,5)
fig, ax = plt.subplots(1)
#ax.set_title('DIFFRACTION', fontsize = 12)
ax.set_xlim([0, 101])
ax.set_ylim([-50, 110])
ax.set_xticks([])
ax.set_yticks([])
plt.axis('off') 

H = 110; h = 18; FS = 12
ax.text(0,H,'wavelength = %2.3e m ' % wL, fontsize = FS, color = [0,0,1])
H = H - h; T = 2*aQx
ax.text(0,H,'Aperature space',fontsize = 12 )
H = H - h    
ax.text(0,H,'  Grid point nQ = %2.0f' % nQ, fontsize = FS, color = [0,0,1])
H = H - h; T = 2*aQx    
ax.text(0,H,'  X width = %2.2e  m' % T, fontsize = FS, color = [0,0,1])
H = H - h; T = 2*aQy  
ax.text(0,H,'  Y width = %2.2e  m' % T, fontsize = FS, color = [0,0,1])
H = H - h  
ax.text(0,H,'Observation ',fontsize = FS )
H = H - h    
ax.text(0,H,'  Grid point nP = %2.0f' % nP, fontsize = FS, color = [0,0,1])
H = H - h    
ax.text(0,H,'  zP = %2.3f   m' % zP, fontsize = FS, color = [0,0,1])
H = H - h    
ax.text(0,H,'First min:  $x_0$ = %2.3f mm' % x0, fontsize = FS, color = [0,0,1])
H = H - h    
ax.text(0,H,'First min:  $y_0$ = %2.3f mm' % y0, fontsize = FS, color = [0,0,1])    
fig.savefig('a4.png')    
    
#%%

fig = plt.figure(figsize =(4, 4))
ax = plt.axes(projection ='3d')
surf = ax.plot_surface(xP, yP, 10*np.log(Irr),
                       cmap = 'hot',
                       edgecolor ='none')
ax.view_init(elev=60, azim=-50, roll=15)
plt.axis('off')
fig.tight_layout()
fig.savefig('a5.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
