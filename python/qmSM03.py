# -*- coding: utf-8 -*-
"""
qmSM02.py    Aug 2024

QUANTUM MECHANICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSM01.pdf


STATISTICAL MECHANICS: MAXWELL - BOLTZMANN DISTRIBUTION
Example: Spectral lines of atomic hydrogen
"""

#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np
import math
from math import factorial 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag
from scipy.integrate import odeint, quad, dblquad, simps

tStart = time.time()

#%%
# Temperature T [K]
T = 50000

# Boltzmann constant
kB = 1.381e-23
# # electron charge
e = 1.602e-19

def fn(c,ER,T):
    g = 2*c**2
    k = -e/(kB*T)
    f = g*exp(k*ER)
    return f

E = zeros(6); ER = zeros(6); nR = zeros(6)
for c in range(6):
    E[c] = -13.6/(c+1)**2
    ER[c] = E[c] - E[0]
    nR[c] = fn(c+1,ER[c],T)
nR = 100*nR/nR[0]

print('Energy levels  [eV]')
print('State    En     ER')
for c in range(6):
    s = c+1; 
    print('  %0.0f' %s + '    %0.2f' %E[c] + '   %0.2f' %ER[c] )
print('  ')
print('Temperature T = %0.0f  K' %T)     
print('Relative populations: nR = 100 for state 1')
print('State     nR ')
for c in range(6):
    s = c+1; print('  %0.0f' %s + '      %0.4f' %nR[c] )






#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)



