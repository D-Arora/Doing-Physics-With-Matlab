# -*- coding: utf-8 -*-
"""
Oct 2026         sim001.py

VIBRATIONS

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/sim001.py

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, odeint, simpson
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, cos, sin 
import time

tStart = time.time()
plt.close('all')

#%% Define parameters
# Grid steps
N = 999
# Max driving force
Fmax = 1
# Frequency range
f1 = 0.1; f2 = 2
# Mass [kg]
m = 0.50666 
# Spring constant [N/m] 
k = 20.0  # Spring constant
# Damping coefficient [kg/s]
b = 2

#%% Calculations
w0 = sqrt(k/m)
f0 = w0/(2*pi)
fD = linspace(f1,f2,N)

wD = 2*pi*fD
A = Fmax / sqrt((k-m*wD**2)**2 + (b*wD)**2)

Find = np.where(A == max(A))
Apeak = A[Find] 
fPeak = fD[Find]

#%% Response 
xP = fD/f0; yP = A
txtX = 'f$_D$ / F$_0$'; txtY = 'A'
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
ax.plot(xP,yP,'b',lw = 2)
fig1.tight_layout()

fig1.savefig('a1.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)


