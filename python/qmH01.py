# -*- coding: utf-8 -*-
"""
qmH01.py    July 2024

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
from sys import exit
pi = math.pi

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
# >>> Principal quantum number: for outputs not computations
n = 1
# >>> ang. mom. quantum number L < n
L = 0
# >>> magnetic quantum number mL <= L             
mL = 0            
# >>> max r value  [m]  enter value in nm x sx
rMax = 2*sx        # r range max [m]

# >>> nuclear charge 
Z = 1

N = 999           # r grid points
num = 100          # number of eigenvalues returned  
rMin = 1e-18       # r range min  [m]


if L > n-1 or abs(mL) > L:
    print('Run aborted: check L < n ans |mL| < L+1')
    exit()
   

#%% SETUP
# >>> max r and min r values for plots: FIG 3
x2 = rMax/sx; x1 = 0             

# Bohr: energies [eV] and radii [nm]
a0 = h**2*eps0/(pi*mu*e**2*Z*sx)
E1 = ( mu*e**4*Z**2/(8*eps0**2*h**2*se) )
 
nq = 25
q = arange(0,nq,1)
#EBohr = zeros(nq); rBohr = zeros(nq)
EBohr = E1/(q+1)**2
rBohr = a0*(q+1)**2

# WARNING: First eignsate found is n = L+! since n > L 
# n > L  -->  adjust principal quantum number for wavefunction plots
nA = n-1
if L > 0:
    nA = nA - L

#%% POTENTIAL ENERGY FUNCTIONS AND PLOTS
r = linspace(rMin,rMax,N)     # radial distance [m]
dx = r[2] - r[1]
# Coulomb potential energy function   [J]
UC = -ke*Z*e**2/r                   
# Orbital kinetic energy function     [J]
UL = L*(L+1)*hbar**2/(2*mu*r**2)
# Effective potential energy function  [J]
U = UC + UL

# FIG 1
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('U  [ eV ]',fontname = 'Times New Roman',fontsize = 14)
ax.set_xlabel('r  [nm]', fontname = 'Times New Roman', fontsize = 14)
ax.set_ylim([-20.1,20.1])
ax.set_xlim([0,1])

ax.set_xlabel('r  [nm]')
ax.set_ylabel('U   [eV]')
ax.set_title('l = %2.0f' % L, style='italic', fontname='Cambria', fontsize = 16)
xP = r/sx
yP = UC/se
ax.plot(xP,yP,'b',lw = 1)
yP = UL/se
ax.plot(xP,yP,'r',lw = 1)
yP = U/se
ax.plot(xP,yP,'k',lw = 2)

fig1.tight_layout()


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
E = real(ev[ev<0]/se)       # Negative eigenvalues [eV]
EB = -E                     # bINDING ENERGIES
lenE = len(E)


#%%  FIG 2: ENERGY PLOTS
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('energy  [ eV ]',fontname = 'Times New Roman',fontsize = 16)
ax.set_xlabel('r  [ nm ]', fontname = 'Times New Roman', fontsize = 16)
ax.set_ylim([-20,2])

xP = [rMin/sx,rMax/sx]
for c in range(lenE):
    yP = [E[c],E[c]]
    ax.plot(xP,yP,'r',lw = 2)
xP = r/sx; yP = U/se
ax.plot(xP,yP,'b',lw = 3)

fig2.tight_layout()


#%%   REDUCED RADIAL WAVEFUNCTION
psi = zeros([N,len(E)]); psi2 = zeros([N,len(E)])
for c in range(len(E)):
    psi[:,c] = ef[:,c]
    psi2[:,c] = psi[:,c]**2
    area = simps(psi2[:,c],r)
    psi[:,c] = psi[:,c]/sqrt(area)

probD = psi**2    # probability density [1/m]

# Eigenstate (ns, L): peak value, normalization, <r>
y = probD[:,nA]
rPeak = r[max(y) == y]/sx

y = psi[:,nA]      # eigenfunction ns
# Probability
fn = y**2
prob = simps(fn,r)

fn = y*r*y
rAvg = simps(fn,r)/sx                 # [nm]
#q = arange(1,lenE,1); 


#%% Prob electron found in region r < a0 expressed as a percentage
q = 0; s = r[0]; c = 0
while s < a0*sx:
    c = c+1
    s = r[c]
    
y1 = psi[0:c,nA]
r1 = r[0:c]
fn = y1**2
prob1 = simps(fn,r1)*100


#%%  CONSOLE OUTPUT
print('  ')
print('grid point N = %2.0f' % N )
print('Z = %2.0f' % Z )
s = rMax/sx;print('rMax = %2.1f nm' % s)
print('ang. mom. quantum no. L = %2.0f' % L)
print('magnetic quantum no. mL = %2.0f' % mL)

if L == 0:
   print(' ')
   print('Z = %2.2f ' %Z + '   Energy [eV]  separation [nm]')
   print('State n    EBohr      Esim        rBohr  ' )
   for q in range(len(E)):
       s = q+1; print(' %0.0f' % s + '           %2.3f' % EBohr[q]
                + '    %2.3f' %EB[q] + '      %2.3f' %rBohr[q])        
print(' ') 

print('Eigenstate: Z = %0.2f  ' %Z +'   n = %0.0f'  %n + '   L = %0.0f' %L 
      +  '   mL = %0.0f' %mL)
s = rMax/sx;print('rMax = %2.1f nm' % s)
print('Percentage probability of electron in region r < a0 = %0.2f' % prob1)
print('EBohr = %0.3f eV' % EBohr[n-1] + '     rBohr = %2.3f nm' % rBohr[n-1])
print('EB    = %0.3f eV' % EB[nA] + '     rPeak = %2.3f nm' % rPeak
      + '    <r> = %2.3f  nm' %rAvg)

 
#%%  FIG 3:  WAVEFUNCTIONS g(R)
# Input principal quantum number, ns and max r value, x2 [nm] for plots 
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,3.2)
fig3, axes = plt.subplots(nrows=1, ncols=2)

C = 0                
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].set_xlabel('r  [ nm ]', fontname = 'Times New Roman', fontsize = 14)
axes[C].set_ylabel('g  ', fontname = 'Times New Roman', fontsize = 14)
axes[C].set_title('n = %2.0f   ' %n  +'  l = %2.0f' % L + '    E = %2.2f eV' %E[nA]
                  , style='italic', 
                  fontname='Cambria', fontsize = 14)
axes[C].set_xlim([x1,x2])
xP = r/sx; yP = psi[:,nA]
if psi[10,0] < 0:
    yP = -yP
yP = yP/ max(yP)    
axes[C].plot(xP,yP, 'blue')
fig3.tight_layout()

C = 1            
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].set_xlabel('r  [ nm ]', fontname = 'Times New Roman', fontsize = 14)
axes[C].set_ylabel('g*g  ', fontname = 'Times New Roman', fontsize = 14)
axes[C].set_xlim([x1,x2])
xP = r/sx; yP = probD[:,nA]
axes[C].plot(xP,yP, 'blue',lw = 2)
xP = [rPeak,rPeak]; yP = [0, max(yP)]
axes[C].plot(xP,yP, 'r',lw =1)
xP = [rAvg,rAvg]; yP = [0, max(yP)]
axes[C].plot(xP,yP, 'm',lw =1)

axes[C].set_title('n = %2.0f   ' % n  +'  l = %2.0f \n' % L +
                  'rPeak = %2.3f ' %rPeak +
                  '   <r> = %2.3f'    %rAvg,
                  style='italic', 
                  fontname='Cambria', fontsize = 12)
fig3.tight_layout()


#%% SAVE FIGURES
#fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)


