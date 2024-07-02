# -*- coding: utf-8 -*-
"""
qm061A.py    June 2024

Ian Cooper 
        matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm061.pdf

QUANTUM MECHANICS     
   SPHERICAL HARMONICS   /   ASSOCIATED LEGENDRE FUNCTIONS
  
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
from scipy.special import clpmn, lpmn
import time

from scipy.special import sph_harm

tStart = time.time()


#%%  INPUTS >>>>>
L = 2             # ORBITAL QUANTUM NUMBER
mL = 1           # MAGNETIC QUANTUM NUMBER
N = 599           # number of grid points
phi = 0           # azimuthal angle  [rad] 
                  # if phi = 0 --> Legendre function


#%%  SETUP
theta = linspace(0,2*pi,N)              # polar angle  [ra]

Y = sph_harm(abs(mL),L,phi,theta).real     # spherical harmonics

probD = Y*Y                           # probability density

#fn = 2*pi*Y*Y*sin(theta)              # check normalization
#prob = simps(fn,theta)                # probability = 1 0r 0
#print('check probability = 1    prob = %2.3f' % prob)


#%% FIG 1:  theta vs Y 
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (3,2)
fig1, ax = plt.subplots(1)

xP = theta/pi; yP = Y
ax.plot(xP,yP,'b',lw=3)    

ax.set_xlim([0,1])

# if prob < 0.9:
#    ax.set_ylim([-1,1])

s = phi/pi
ax.set_title('L = %2.0f' % L + '    m$_L$ = %2.0f' %mL +
            '   $\phi$ / $\pi$ = %2.2f' % s)
ax.set_ylabel('Y',color= 'b')    
ax.set_xlabel(r'$\theta$ / $\pi$ ') 
ax.grid()       
fig1.tight_layout()       


#%% FIG 2:  probability density
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (3,2)
fig2, ax = plt.subplots(1)

xP = theta/pi; yP = probD


if mL < 0:
     Y = sqrt(2) * (-1)**mL * Y.imag
elif mL > 0:
     Y = sqrt(2) * (-1)**mL * Y.real
yP = Y
ax.plot(xP,yP,'b',lw=3) 

#ax.set_ylim([-1,1])
ax.set_xlim([0,1])
s = phi/pi
ax.set_title('L = %2.0f' % L + '    m$_L$ = %2.0f' %mL +
             '   $\phi$ / $\pi$ = %2.2f' % s)
ax.set_ylabel('Y*Y',color= 'b')    
ax.set_xlabel(r'$\theta$ / $\pi$ ') 
ax.grid()       
fig2.tight_layout()    


#%% FIG 3  probability density   polar plot
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (3,3)
fig3, ax = plt.subplots(nrows=1, ncols=1,subplot_kw={'projection': 'polar'})


ax.plot(theta,probD,'b',lw=3)
#ax.plot(theta,Y,'r',lw=3)


ax.set_rticks([]) 

ax.grid(True)
ax.set_xticks(arange(0,2*pi, pi/6))
ax.set_title('L = %2.0f' % L + '    m$_L$ = %2.0f' %mL +
             '   $\phi$ / $\pi$ = %2.2f' % s)
fig3.tight_layout() 

#%% SAVE FIGURES
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)


