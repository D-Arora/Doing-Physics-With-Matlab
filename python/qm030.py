# -*- coding: utf-8 -*-
"""
qm030.py            May 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

QUANTUM MECHANICS
   Scattering from a finite setp potential
   Numerical solutions
   
   https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm030.pdf

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
from scipy.integrate import odeint, solve_ivp



#%% FUNCTIONS
def lorenz(x, state):    
    u, v = state
    P = E0
    
    if x > 0:
       P = E0 - U0
    
    du = v 
    dv = C*P*u
    return [du, dv]


#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [ J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
me = 9.10938356e-31            # Electron mass [kg]

xMin = -0.6e-9; xMax = 0.6e-9     # X domain  [m]
E0 = 80                   # beam particle energy  [eV]
U0 = 50             # height of potential step [eV]
N = 8299                           # Grid points

            
#%% COMPUTATIONS
x = np.linspace(xMin,xMax,N)
E = e*E0                          # beam particle energy  [J]                                 
U = e*U0                          # height of potential step [J]
w = E/hbar                        # angular frequency [rad/s]
T = 2*pi/w                        # oscillation period [s]
k1 = sqrt(2*me*E)/hbar            # propagation constant x < 0    [1/m]
k2 = cmath.sqrt(2*me*(E-U))/hbar       # propagation constant x > 0    [1/m]
L1 = 2*pi/k1                      # wavelength x < 0    [m]
L2 = 2*pi/k2                      # wavelength x > 0    [m]

C = -2*me*e/hbar**2                 # SE constant        

# SOLVE SCHRODINGER EQUATION
xSpan = np.flip(x)

u1 = 0; u2 = k1
u0 = [u1,u2]
sol = odeint(lorenz, u0, xSpan,  tfirst=True)
psiR = sol[:,0]

u1 = 1; u2 = 0
u0 = [u1,u2]
sol = odeint(lorenz, u0, xSpan, tfirst=True)
psiI = sol[:,0]

psiR = np.flip(psiR)/np.amax([psiI,psiR])
psiI = np.flip(psiI)/np.amax([psiI,psiR])

probD = psiR**2 + psiI**2


#%% GRAPHICS
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6,4)

fig1, axes = plt.subplots(nrows=2, ncols=1)
xP = x*1e9

R = 0
axes[R].set_ylabel('$\psi$',color = 'black',fontsize = 14)
axes[R].xaxis.grid()
axes[R].yaxis.grid()
axes[R].set_xlim(xMin*1e9,xMax*1e9)
axes[R].set_xticks(np.arange(-0.6,0.7,0.2))
yP = psiR
axes[R].plot(xP, yP,'b', lw = 2)
yP = psiI
axes[R].plot(xP, yP,'r', lw = 2)

R = 1
axes[R].set_xlabel('x  [ nm ] ',color = 'black',fontsize = 10)
axes[R].set_ylabel('|$\psi|^2$',color = 'black',fontsize = 14)
axes[R].xaxis.grid()
axes[R].yaxis.grid()
axes[R].set_xlim(xMin*1e9,xMax*1e9)
yP = probD
axes[R].plot(xP, yP,'b', lw = 2)
axes[R].fill_between(xP, probD,color = [1,0,1],alpha=0.2)

fig1.tight_layout()

fig1.savefig('a3.png')




       