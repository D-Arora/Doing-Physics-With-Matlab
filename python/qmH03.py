# -*- coding: utf-8 -*-
"""
qm050.py    June 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm050.pdf

QUANTUM MECHANICS     HCl molecule
   Solution of the time independent Schrodinger equation
   for an electron bound in a potential well by finding the eigenvalues
   and eigenvectors
   Parabolic potential well and Morse potential


https://physicspython.wordpress.com/tag/hydrogen-atom/

"""

# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt
import scipy.special 
import scipy.special as sps
from scipy import *
import scipy.linalg as la
import math
import cmath
import os
import matplotlib.pyplot as plt
from pylab import *
from math import *
from numpy import linspace,matrix,array
from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve
from scipy.integrate import odeint, solve_ivp, simps
from scipy.sparse.linalg import eigsh, eigs #Solves the Eigenvalue problem
from scipy.sparse import diags #Allows us to construct our matrices
from matplotlib.animation import FuncAnimation, PillowWriter 
import time
from scipy.special import sph_harm

tStart = time.time()


#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
cL = 2.99792458e8              # speed of light
eps0 = 8.8541878188e-12        # episilon-zero  [F/m]
me = 9.109e-31
mp = 1.673e-27
mu = me*mp/(me+mp)             # reduced mass HCl molecule [kg]
se = e                         # Energy scaling factor   J <---> ev 
sx = 1e-9                      # Length scaling factor   m <---> nm
ke = 1/(4*pi*eps0)             # Coulomb constant

#%%  INPUTS
n = 4       # principle quantum number  n > L
L =  3      # orbital ang. mom. quantum number L == l: L < n
mL = 3     # magnetic quantum number mL <= L

rMax = 6*sx       # r range max [m]

Z = 1               # nuclear charge
N = 255             # r grid points
num = 100           # number of eigenvalues returned  
rMin = 1e-18        # r range min  [m]


#%% SETUP
#n > L  -->  adjust principal quantum number for wavefunction plots
ns = n 
if L > 0:
   n = n - L

r = linspace(rMin,rMax,N)     # radial distance [m]
dx = r[2] - r[1]
# Coulomb potential energy function   [J]
UC = -ke*Z*e**2/r                   
# Orbital kinetic energy function     [J]
UL = L*(L+1)*hbar**2/(2*mu*r**2)
# Effective potential energy function  [J]
U = UC + UL

#%% SPHERICAL HARMONICS

theta = linspace(0,2*pi,N)
Y = sph_harm(abs(mL), L, 0, theta)
YR = Y.real; YI = Y.imag

if mL < 0:
      Y = sqrt(2) * (-1)**mL * Y.imag
elif mL > 0:
      Y = sqrt(2) * (-1)**mL * Y.real
          

#%% EIGENVALUES AND EIGENFUNCTIONS
Cse = -hbar**2/(2*mu)               # Schrodinger constant
UM = diag(U)                        # potential energy matrix
               
# AM (second derivative), KM (kinetic energy), HM (Hamiltonian) matrices
off = ones(N-1)
AM = (-2*np.eye(N) + np.diag(off,1) + np.diag(off,-1))/(dx**2)
KM = Cse*AM
HM = KM + UM

# Eigenvalues [J] and eigenfunctions (eigenvectors)
ev, ef = eigsh(HM, which="SM", k = num)

# OUTPUT: negative eigenvalues and normalize eigenfunctions
E = real(ev[ev<0]/se)                 # negative eigenvalues [eV]
EB = E
lenE = len(E)
print('  ')
print('Numerical',EB)
print('  ')


#%%   REDUCED RADIAL WAVEFUNCTION
psi = zeros([N,len(E)]); psi2 = zeros([N,len(E)])
for c in range(len(E)):
    psi[:,c] = ef[:,c]
    psi2[:,c] = psi[:,c]**2
    area = simps(psi2[:,c],r)
    psi[:,c] = psi[:,c]/sqrt(area)

g = psi[:,n-1]
probD = psi**2    # probability density [1/m]


#%% FIG 1: Reduced radial function and spherical harmonics
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6.2,3.2)
fig1, axes = plt.subplots(nrows=1, ncols=2)

C = 0                
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].set_xlabel('r  [ nm ]', fontname = 'Times New Roman', fontsize = 14)
axes[C].set_ylabel('g(r)*g(r)  ', fontname = 'Times New Roman', fontsize = 14)
axes[C].set_title('n = %2.0f   ' %ns  + '  l = %2.0f' % L 
             + '    ml = %2.0f' % mL 
             + '    E = %2.2f eV' %E[n-1],
             style='italic', 
             fontname='Cambria', fontsize = 14)
#axes[C].set_xlim([x1,x2])
xP = r/sx; yP = g*g
if psi[10,0] < 0:
    yP = -yP
yP = yP/ max(yP)    
axes[C].plot(xP,yP, 'blue')
fig1.tight_layout()

C = 1            
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].set_xlabel(r'$\theta / \pi$', fontsize = 12)
axes[C].set_ylabel(r'$Y(\theta)*Y(\theta)$  ', fontsize = 12)
xP = theta/pi; yP = Y*Y
axes[C].plot(xP,yP, 'blue',lw = 2)
fig1.tight_layout()

#%% FIG 2: PROBABILITY DENSITY  pd
x = zeros([N,N])
z = zeros([N,N])
pd = zeros([N,N])

for c1 in range(N):
    for c2 in range(N):
        x[c1,c2] = r[c2]*sin(theta[c1])
        z[c1,c2] = r[c2]*cos(theta[c1])
        pd[c1,c2] = g[c2]**2 * Y[c1]*Y[c1]        

pd = pd/amax(pd)
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3.2)
rcParams["patch.edgecolor"]
fig2, ax = plt.subplots(nrows=1, ncols=1)

plt.pcolormesh(x/sx,z/sx,pd**0.5,cmap='jet',edgecolors ='none',)

ax.set_xlabel('x  [ nm ]', fontname = 'Times New Roman', fontsize = 14)
ax.set_ylabel('z  [ nm ]  ', fontname = 'Times New Roman', fontsize = 14)
ax.set_title('n = %2.0f   ' %ns  + '  l = %2.0f' % L 
             + '    ml = %2.0f' % mL 
             + '    E = %2.2f eV' %E[n-1],
             style='italic', 
             fontname='Cambria', fontsize = 14)

ax.set_aspect('equal', adjustable='box')
col = [0,0,1]
ax.set_facecolor((0,0,0.5))
fig2.tight_layout()


#%%  FIG 3: SURFACE PLOT
plt.rcParams["figure.figsize"] = (5,3.2)
fig3 = plt.figure(figsize=plt.figaspect(1.))
ax = fig3.add_subplot(111, projection='3d')
ax.plot_surface(x/sx, z/sx, pd**0.5, cmap = 'jet' )
ax.set_axis_off()


#%% SAVE FIGURES
fig1.savefig('a1.png')
# fig2.savefig('a2.png')
# fig3.savefig('a3.png')

fig2.savefig('a2.png',bbox_inches='tight',dpi=800,
             pad_inches=0.1, transparent=True)
fig3.savefig('a3.png',bbox_inches='tight',dpi=800,
             pad_inches=0.1, transparent=True)


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)


