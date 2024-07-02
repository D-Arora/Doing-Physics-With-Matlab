# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:13:47 2024

@author: Owner
"""

"""
qm061C.py    June 2024

Ian Cooper 
        matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm061.pdf

QUANTUM MECHANICS     
   SPHERICAL HARMONICS   /   ASSOCIATED LEGENDRE FUNCTIONS
   
Origin source of Code
https://scipython.com/blog/visualizing-the-real-forms-of-the-spherical-harmonics/

https://people.csail.mit.edu/sparis/sh/index.php?img=64
   
"""

# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, radians
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
from scipy.special import sph_harm
import time

tStart = time.time()


#%% INPUTS
# oribtal quantum nuumber L (L == l)
L = 1
# magnetic quantum number mL   (mL == ml)
mL = -1
# Grid point N
N = 299

 
#%% SETUP
# Polar angle  [rad]  / Azimuthal angle  [ra]
theta = linspace(0,pi, N)
phi   = linspace(0,2*pi, N)
# [D] meshgrid
THETA, PHI =  np.meshgrid(theta, phi) 
# Calculate the Cartesian coordinates of each point in the mesh
xyz = np.array([sin(THETA)*sin(PHI), sin(THETA)*cos(PHI), cos(THETA)])
# Calculate spherical harmonic
Y = sph_harm(abs(mL), L, PHI, THETA)
YR = Y.real; YI = Y.imag

if mL < 0:
     Y = sqrt(2) * (-1)**mL * Y.imag
elif mL > 0:
     Y = sqrt(2) * (-1)**mL * Y.real

Yx, Yy, Yz = abs(Y) * xyz


plt.rcParams["figure.figsize"] = (2,2)
fig3 = plt.figure()
fig3.subplots_adjust(top = 0.8, bottom = 0.0, left = 0.0,
                    right = 1, hspace = 0.20,wspace=0.20)
 
# Creating our empty subplots
ax = fig3.add_subplot(1, 1, 1, projection='3d')
 


ax.plot_surface(Yx, Yy, Yz, cmap='jet')
ax.set_title('l = %2.0f' %L + '   $m_l = %2.0f$' %mL)
L2 =0.2; L1 = -L2
ax.plot([L1,L2], [0,0], [0,0], c='0.4', lw=1, zorder=10)
ax.plot([L1,L2], [0,0], [0,0], c='0.4', lw=1, zorder=10)
ax.plot([0,0], [L1,L2], [0,0], c='0.4', lw=1, zorder=10)
ax.plot([0,0], [0,0], [L1,L2], c='0.4', lw=1, zorder=10)
ax.set_aspect('equal')

ax.axis('off')

#fig3.savefig('a1.png',bbox_inches='tight',dpi=800,
#             pad_inches=0.1, transparent=True)