# -*- coding: utf-8 -*-

"""
qm061.py    June 2024

Ian Cooper 
        matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm061.pdf

QUANTUM MECHANICS     ANGULAR MOMENTUM:  AZIMUTHAL WAVEFUNCTION
  
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
import time

tStart = time.time()


#%%
phi = linspace(0,2*pi,599)         # arimuthal angle  [rad]

# r real part / im imaginary part / pd probability density

def graph(R,C,mL):
    r = abs(cos(mL*phi)); im = abs(sin(mL*phi)); pd = r**2 + im**2
    ax[R,C].plot(phi, r,'b')
    ax[R,C].plot(phi, im,'k')
    ax[R,C].plot(phi, pd,'r',lw = 4)
    ax[R,C].set_rmax(1)
    ax[R,C].set_rticks([])  # Less radial ticks
    ax[R,C].grid(True)
    ax[R,C].set_xticks(arange(0,2*pi, pi/6))
    ax[R,C].set_title('m$_L$ = %2.0f' % mL)
    

#%%
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7.2)
fig1, ax = plt.subplots(nrows=3, ncols=2,subplot_kw={'projection': 'polar'})
#fig1, ax = plt.subplots(subplot_kw={'projection': 'polar'})
fig1.subplots_adjust(top = 0.95, bottom = 0.07, left = 0.120,\
                    right = 0.98, hspace = 0.60,wspace=0.40)

R = 0; C = 0; mL = 0 
graph(R,C,mL)
R = 0; C = 1; mL = 1  
graph(R,C,mL)
R = 1; C = 0; mL = 2  
graph(R,C,mL)
R = 1; C = 1; mL = 3  
graph(R,C,mL)
R = 2; C = 0; mL = 4 
graph(R,C,mL)
R = 2; C = 1; mL = 5  
graph(R,C,mL)

fig1.tight_layout() 
#fig1.savefig('a1.png')




#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)





