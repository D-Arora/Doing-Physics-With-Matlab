# -*- coding: utf-8 -*-
"""
op005.py               March 2024
 
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
nQ = 599              # aperture points for Q  must be ODD

wL = 632.8e-9        # wavelength [m]

# Aperautre space: full-width [m]
aQx = 1e-4; aQy = 2e-4           
 
#Observation spacehalf:  zP   [m] 
zPmin = 1e-5; zPmax = 100e-5;       


         

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
EP = np.zeros(nP)+np.zeros(nP)*1j         

zP = linspace(zPmin,zPmax,nP)     



# Simpson [2D] coefficients
S = np.ones(nQ)
R = np.arange(1,nQ,2);   S[R] = 4;
R = np.arange(2,nQ-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy

# % COMPUATION OF DIFFRACTION INTEGRAL FOR ELECTRIC FIELD % IRRADIANCE -----

for c1 in range(nP):
         rPQ = np.sqrt(xQ**2 + yQ**2 + zP[c1]**2)
         rPQ3 = rPQ*rPQ*rPQ
         kk = ik * rPQ
         MP1 = exp(kk)
         MP1 = MP1 / rPQ3
         MP2 = zP[c1] * (ik * rPQ - unit)
         MP = MP1 * MP2
         EP[c1] = sum(sum(MP*S))

Irr = np.real(EP*np.conj(EP))
Irr = Irr/amax(amax(Irr))


#%%  GRAPHICS --------------------------------------------------------------
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3.9,3)


fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(zP*1e3,Irr,'r',lw = 2)
ax.set_xlabel('$z_P$  [mm]', fontsize=12)
ax.set_ylabel('I  [a.u.]', fontsize=12)
#ax.set_xlim([-20,20])
ax.xaxis.grid()
ax.yaxis.grid()
#ax.legend()
fig.tight_layout()
fig.savefig('a1.png')

fig, ax = plt.subplots(nrows=1, ncols=1)
plt.yscale('log');plt.xscale('log')
ax.plot(zP,Irr,'r',lw = 2)
ax.set_xlabel('$z_P$   [mm]', fontsize=12)
ax.set_ylabel('I  [a.u.]', fontsize=12)
#ax.set_xlim([-20,20])
ax.xaxis.grid()
ax.yaxis.grid()
ax.axis('equal')
#ax.legend()
fig.tight_layout()
fig.savefig('a2.png')





#%%   NUMERICAL SUMMARY
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3.3,5)
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


  
fig.savefig('a4.png')    
    




#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
