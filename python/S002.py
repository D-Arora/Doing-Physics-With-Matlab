# -*- coding: utf-8 -*-
"""
S001.py               27 May 2025
 
 RAY OPTICS
 Matrix methods in paraxial optics
 Paraxial approximation
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/S002.pdf

"""

import numpy as np
from numpy import pi, arange, exp, sin, linspace, zeros, amax, sqrt, tan
import matplotlib.pyplot as plt
import time
from scipy import special
from scipy.signal import find_peaks
from scipy.special import jv
from mpl_toolkits.mplot3d import Axes3D
tStart = time.time()

plt.close('all')


#%% PARAXIAL APPROXIMATION
# Angle alpha == A
A1 = 0; A2 = 20; N = 999
A = linspace(A1,A2,N)        # [deg]
AR = A*pi/180                # [rad]

m = tan(AR)                  # slope 
E = 100*(m-AR)               # percentaage error
AE = A[E>1][0]               # Angle A when E = 1%

#%% GRAPICS
fig1, ax = plt.subplots(nrows=1, ncols=1)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)

ax.plot(A,m,'b',lw = 2,label = r'tan$\alpha$')
ax.plot(A,AR,'r',lw = 2, label = r'$\alpha$')

ax.set_xlabel(r'angle $\alpha$  [deg]', fontsize=12)
ax.set_ylabel('slope  m', fontsize=12)
ax.set_title('paraxial approx for slope  AE = %0.1f deg' %AE,fontsize = 10)
ax.legend(fontsize = 10)
ax.grid()

ax1 = ax.twinx()
ax1.plot(A,E,'k',lw = 1,)
ax1.plot([AE,AE], [0,1],'m',lw = 1)
ax1.set_ylabel('E %', fontsize=12)
fig1.tight_layout()

fig1.savefig('a1.png')


