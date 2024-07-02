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
from scipy.special import sph_harm
tStart = time.time()


#%%
N = 299
theta = linspace(0,pi, N)
phi   = 0

#%%
# r real part / im imaginary part / pd probability density

def graph(R,C,L,mL):
    Y = sph_harm(abs(mL), L, mL, theta)

    if mL < 0:
         Y = np.sqrt(2) * (-1)**mL * Y.imag
    elif mL > 0:
         Y = np.sqrt(2) * (-1)**mL * Y.real

    r = conj(Y)*Y.real
    
    ax[R,C].plot(theta+pi/2, r,'b',lw = 2)
    ax[R,C].plot(-theta+pi/2, r,'r',lw = 2)
    #ax[R,C].set_rmax(1)
    ax[R,C].set_rticks([])  
    ax[R,C].set_theta_offset(pi/2)
    ax[R,C].grid(True)
    ax[R,C].set_xticks(arange(0,2*pi, pi/6))
    ax[R,C].set_title('L = %2.0f' % L + '   m$_L$ = %2.0f' % mL)



#%%
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7.2)
fig1, ax = plt.subplots(nrows=3, ncols=2,subplot_kw={'projection': 'polar'})

LA = array([0.,1.,1.,1.,2.,2.])
MA = array([0.,-1.,0.,1.,-2.,-1.])

# LA = array([2.,2.,2.,3.,3.,3.])
# MA = array([0.,1.,2.,-3.,-2.,-1.])

LA = array([3.,3.,3.,3.,4.,4.])
MA = array([0.,1.,2.,3.,0.,1.])

  
R = 0; C = 0; L = LA[0];mL= MA[0] 
graph(R,C,L,mL)
R = 0; C = 1; L = LA[1];mL= MA[1]
graph(R,C,L,mL)
R = 1; C = 0; L = LA[2];mL= MA[2]
graph(R,C,L,mL)
R = 1; C = 1; L = LA[3];mL= MA[3]
graph(R,C,L,mL)
R = 2; C = 0; L = LA[4];mL= MA[4]
graph(R,C,L,mL)
R = 2; C = 1; L = LA[5];mL= MA[5]
graph(R,C,L,mL)

fig1.tight_layout() 






#%% SAVE FIGURES
fig1.savefig('a1.png')




#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)





