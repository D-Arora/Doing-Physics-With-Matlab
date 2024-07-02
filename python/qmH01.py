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
a0 = 5.292e-11                  # Bohr radius  [m]
ER = hbar**2/(2*mu*a0**2)/se

#%%  INPUTS
L = 3              # orbital ang. mom. quantum number L == l: L < n
mL = 1             # magnetic quantum number

rMax = 5*sx        # r range max [m]

ns = 4             #  FIG 3: Principal quantum number for wavefunction plots
x2 = 3             #  FIG 3: max r values for wavefunction plots  [nm]

Z = 1              # nuclear charge
N = 1201           # r grid points
num = 100          # number of eigenvalues returned  
rMin = 1e-18       # r range min  [m]
x1 = 0             # min r value for  wavefunction plots  [nm]

#%% SETUP
# n > L  -->  adjust principal quantum number for wavefunction plots
nsp = ns
if L > 0:
    ns = ns - L

r = linspace(rMin,rMax,N)     # radial distance [m]
dx = r[2] - r[1]
# Coulomb potential energy function   [J]
UC = -ke*Z*e**2/r                   
# Orbital kinetic energy function     [J]
UL = L*(L+1)*hbar**2/(2*mu*r**2)
# Effective potential energy function  [J]
U = UC + UL

# Theoretical energy eigenvalues [eV] and Bohr radii  [nm]
ne = arange(1,20,1)
ET = ER/ne**2
rB = (hbar**2/(sx*ke*mu*e**2))*ne**2
a0 = rB[0]      # Bohr radius

#%% FIG 1: Potential energy functions
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('U  [ eV ]',fontname = 'Times New Roman',fontsize = 14)
ax.set_xlabel('r  [nm]', fontname = 'Times New Roman', fontsize = 14)
ax.set_ylim([-10.1,2.1])
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
E = real(ev[ev<0]/se)                 # negative eigenvalues [eV]
EB = E
lenE = len(E)
print("Energies E  [eV]")
print('Theoretical',ET)
print('  ')
print('Numerical',EB)
print('  ')
print('Bohr radii  [nm]')
print(rB)

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
y = probD[:,ns-1]
rPeak = r[max(y) == y]/sx

y = psi[:,ns-1]      # eigenfunction ns
# Probability
fn = y**2
prob = simps(fn,r)
# Position and its uncertainty
fn = y*r*y
rAvg = simps(fn,r)/sx                 # [nm]


#%%  FIG 3:  WAVEFUNCTIONS g(R)
# Input principal quantum number, ns and max r value, x2 [nm] for plots 
# ns = 2
# x2 = 2
# x1 = 0
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,3.2)
fig3, axes = plt.subplots(nrows=1, ncols=2)

C = 0                
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].set_xlabel('r  [ nm ]', fontname = 'Times New Roman', fontsize = 14)
axes[C].set_ylabel('g  ', fontname = 'Times New Roman', fontsize = 14)
axes[C].set_title('n = %2.0f   ' %nsp  +'  l = %2.0f' % L + '    E = %2.2f eV' %E[ns-1]
                  , style='italic', 
                  fontname='Cambria', fontsize = 14)
axes[C].set_xlim([x1,x2])
xP = r/sx; yP = psi[:,ns-1]
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
xP = r/sx; yP = probD[:,ns-1]
axes[C].plot(xP,yP, 'blue',lw = 2)
xP = [rPeak,rPeak]; yP = [0, max(yP)]
axes[C].plot(xP,yP, 'r',lw =1)
xP = [rAvg,rAvg]; yP = [0, max(yP)]
axes[C].plot(xP,yP, 'm',lw =1)

axes[C].set_title('n = %2.0f   ' % nsp  +'  l = %2.0f \n' % L +
                  'rPeak = %2.3f ' %rPeak +
                  '   <r> = %2.3f'    %rAvg,
                  style='italic', 
                  fontname='Cambria', fontsize = 12)


fig3.tight_layout()

# #%% CONSOLE OUTPUT
# print('  ')
# print('grid point N = %2.0f' % N + '   eigenvalues returned M = %2.0f' % num)
# s1 = xMin/sx; s2 = xMax/sx; print('xMin = %2.2f nm' % s1 + '   xMax = %2.2f nm' % s2)
# s = x0/sx; print('Equilibrium bond length x0 = %2.3f nm' %s)  
# s = U0/se; print('well depth U0 = %2.3f ev' %s) 
# print('spring constant k = %2.3f N/m' % k) 
# print(' ')
# print('Energies [eV]')
# print('State n    Ew       ET      dE       E       EB')
# for q in range(len(E)):
#     print(' %0.0f' % q + '        %2.3f' % Ew[q] + '    %2.3f' %ET[q] 
#           + '   %2.3f' % dE[q] + '   %2.3f' %E[q] +  '   %2.3f'  %EB[q])
# print(' ')  




# #%% SCATTER PLOT FOR PROBILITY
# q = n
# y = probD[:,q]/amax(probD)

# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (6,1)
# fig11, ax = plt.subplots(1)
# fig11.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.18,\
#                     right = 0.80, hspace = 0.20,wspace=0.2)

# x1 = xMin/sx; x2 = xMax/sx
# ax.set_xlim([x1,x2])
# ax.set_ylim([0,1])
# ax.set_axis_off()

# random.seed()
# num = 20000
# for c in range(num):
#     xP = x1+x2*random.random()        # [mm]
#     yP = random.random()
#     qR = random.random()
#     qP = y[x/sx>xP]
#     q0 = qP[0] 
#     if qR < q0:
#        ax.plot(xP,yP,'bo', ms = 1)
# fig11.savefig('a11.png')
#%% SAVE FIGURES
#fig1.savefig('a1.png')
#fig2.savefig('a2.png')
fig3.savefig('a3.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)


