# -*- coding: utf-8 -*-
"""
qmCB1.py    July 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmCB.pdf

QUANTUM MECHANICS     
     Double well potential: COVALNET BONDING
     
TAKES A FEW MINUTES TO RUN    ~ 8 minutes
Reduce N to reduce run time     
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
from scipy.interpolate import make_interp_spline
import time

tStart = time.time()


#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
cL = 2.99792458e8              # speed of light
me = 9.10938356e-31           # electron mass [kg]
eps0 = 8.8541878188e-12        # episilon-zero  [F/m]
se = e                         # Energy scaling factor   J <---> ev 
sx = 1e-9                      # Length scaling factor   m <---> nm
ke = 1/(4*pi*eps0)             # Coulomb constant


#%%  INPUTS
xw = 0.5                   # width of of left well [nm]
xb = 0.5                   # width of boundaries [nm] 
U0 = -50                   # depth of well

N = 2599                   # x grid points
num = 50                   # number of eigenvalues returned  


#%% SETUP: potential energy function
ns = 30
xs = linspace(0.001,0.05,ns)
EA = zeros(ns);ES = zeros(ns)
for q in range(ns):
    L = xb + xw + xs[q] + xw + xb       # x range  [nm]
    x = linspace(-L/2, L/2, N)          # x grid   [nm]

# Double well
    U = zeros(N)                       # potential energy  [eV]
    X1 = xs[q]/2; X2 = X1+xw
    U[x>X1] = U0 
    U[x>X2] = 0
    U[x<-X1] = U0
    U[x<-X2] = 0                       # potential energy [J}]
    Ue = U*se
    xe = x*sx
    dx = xe[2] - xe[1] 

# EIGENVALUES AND EIGENFUNCTIONS
    Cse = -hbar**2/(2*me)                   
    UM = diag(Ue)  
               
# AM (second derivative), KM (kinetic energy), HM (Hamiltonian) matrices
    off = ones(N-1)
    AM = (-2*np.eye(N) + np.diag(off,1) + np.diag(off,-1))/(dx**2)
    KM = Cse*AM
    HM = KM + UM

# Eigenvalues [J] and eigenfunctions (eigenvectors)
    ev, ef = eigsh(HM, which="SM", k = num)
    
# OUTPUT: negative eigenvalues and normalize eigenfunctions
    E = ev[ev<0]/se 
    ES[q] = E[0]
    EA[q] = E[1]                # negative eigenvalues [eV]
    

#%% RESULTS for all separation distances
E0 = -48.77690
psiS = ef[:,0]
psiA = ef[:,1]


#%%  FIG 1: ENERGY PLOTS
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(1)

xP = xs; yP = ES-E0; zP = EA-E0

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('energy  [ eV ]',color= 'black')
ax.set_xlabel('separation distance xs  [nm]',color = 'black')
# # #ax.set_ylim([-5,0])
ax.plot(xP,yP,'bo',ms = 3,label = 'symmetric')
ax.plot(xP,zP,'ro',ms = 3,label='antisymmetric')
ax.legend()
fig1.tight_layout()


#%%  NUCLEAR REPULSION
A = 1.5e-3
# xR = linspace(0.01,0.12,199)
ER = A/xs
Y = ER+(ES-E0)
X_Y_Spline = make_interp_spline(xs, Y)
X_ = np.linspace(xs.min(), xs.max(), 500)
Y_ = X_Y_Spline(X_)
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('energy  [ eV ]',color= 'black')
ax.set_xlabel('separation distance xs  [nm]',color = 'black')
# # #ax.set_ylim([-5,0])
ax.plot(xP,yP,'bo',ms = 3,label = 'symmetric')
ax.plot(xP,zP,'ro',ms = 3,label='antisymmetric')
#ax.plot(xs,ER,'bo',ms = 1)
xRP = linspace(0.002,0.05,199); yRP = A/xRP
ax.plot(xRP,yRP,'k',lw = 1,label = 'nuclear repulsion')

ax.plot(xs,Y,'mo',ms = 3, label = 'total energy')

ax.plot(X_, Y_,'m',lw = 2)

ax.legend()
fig2.tight_layout()



#%%  SAVE FIGURES
# fig1.savefig('a1.png')
# fig2.savefig('a2.png')
# fig3.savefig('a3.png')
# fig4.savefig('a4.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)

