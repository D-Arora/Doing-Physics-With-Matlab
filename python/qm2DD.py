# -*- coding: utf-8 -*-
"""
qm2DD.py            July 2024

Ian Cooper
      email: matlabvisualphysics@gmail.com

Documentation
     https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm2Dfdtd.pdf

QUANTUM MECHANICS
   Finite Difference Time Development Method
   [2D] Schrodinger Equation
   Propagation of a [2D] Gausssian Pulse
   Arbitrary units are used
   DOUBLE SLIT

Code - modified version(Kevin Berwick)
        'Computational Physics' by Giordano Nakanishi

"""

#%% LIBRARIES

import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax 
from numpy.linalg import eig
from scipy.integrate import odeint, quad, dblquad, simps
from scipy import pi, sqrt
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import axes3d
import sympy as sym
import sys

tStart = time.time()


#%% SETUP XY GRID
# Number of time Steps [400]
nT = 600
# Mesh points X and y  [180] and XY grid
N = 208
G = linspace(0,1,N)
xG, yG = np.meshgrid(G,G)
dx = G[2] - G[1]
# Constant in S.E. and time step
f = 0.2             # f < 0.2
dt = 2*dx**2*f


#%% Potential energy function
U0 = 25000
N2 = int(N/2)
U = zeros([N,N])
U[0:N2-10, N2-5:N2+5]  = U0
U[N2+10:N, N2-5:N2+5]  = U0
U[N2-5:N2+5, N2-5:N2+5] = U0
sf = 0.1    


#%%
# [2D] GAUSSIAN PULSE (WAVE PACKET) 
# Initial centre of pulse
x0 = 0.20;  y0 = 0.5
# Wavenumber
k0 = 100
# Initial amplitude of pulse 
A = 10
# Pulse width: sigma squared
s = 10e-3
# Envelope
psiE = A*exp(-(xG-x0)**2/s)*exp(-(yG-y0)**2/s)
# Plane wave propagation in +X direction
psiP = exp(1j*k0*xG)
# Wavefunction
psi1 = psiE*psiP
# Probability Density  
prob1 = np.conj(psi1)*psi1
# Extract Real and Imaginary parts
R1 = np.real(psi1);  I1 = np.imag(psi1)


#%% FUNCTIONS
def fI(I1,R1):
  I2 = zeros([N,N])
  for x in range(2,N-1,1):
    for y in range(2,N-1,1):
       I2[x,y] = I1[x,y] + f*(R1[x+1,y] - 2*R1[x,y] + R1[x-1,y]
                         +    R1[x,y+1]  -2*R1[x,y] + R1[x,y-1]) - dt*U[x,y]*R1[x,y]  
  return I2 

def fR(I1,R1):
  R2 = zeros([N,N])
  for x in range(2,N-1,1):
    for y in range(2,N-1,1):
       R2[x,y] = R1[x,y] - f*(I1[x+1,y] - 2*I1[x,y] + I1[x-1,y]
                         +    I1[x,y+1] - 2*I1[x,y] + I1[x,y-1]) + dt*U[x,y]*I1[x,y]  
  return R2 
 

#%%   SOLVE [2D] Schrodinger equation 
for c in range(1,nT,1):
# Update real part of wavefunction    
  R2 = fR(I1,R1)
  R1 = R2

# Update imaginary part of wavefunction
  I2 = fI(I1,R1)

# Probability Density Function  
  probD = R2**2 + I1**2
  I1 = I2


#%% FIG 1: PROBABILITY DENSITY: PCOLOR plot at end of simulation
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6,2.3)
fig1, ax = plt.subplots(1)

ax.pcolormesh(xG,yG,probD**0.1)
ax.plot(x0,y0,'ro')
ax.plot([0.5,0.5],[0.3,0.7],'r',lw =3)
ax.set_ylim([0.3,0.7])
ax.set_aspect('equal', adjustable='box')
plt.text(0.05,0.65, 'nT',color = 'yellow', fontsize=12 )
plt.text(0.05,0.60, nT,color = 'yellow',fontsize = 12 )
plt.text(0.05,0.38, '$k_0$',color = 'yellow', fontsize=12 )
plt.text(0.05,0.33, k0,color = 'yellow',fontsize = 12 )
plt.text(0.85,0.65, '$U_0$',color = 'yellow', fontsize=12 )
plt.text(0.85,0.60, U0,color = 'yellow', fontsize=12 )
ax.axis('off')
fig1.tight_layout()


#%%  FIG 3: POTENTIAL ENERGY [3D] plot
fig3, ax = plt.subplots(1)
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3,3)
#surf = ax.pcolormesh(xG,yG,U)
surf = ax.pcolormesh(xG,yG,U)
ax.set_xticks([0,0.5,1])
ax.set_yticks([0,0.5,1])
ax.set_xlabel('x'); ax.set_ylabel('y')
ax.set_aspect('equal', adjustable='box')
fig3.tight_layout()
 

#%% SAVE FIGURES
# fig1.savefig('a1.png')
# fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)   
