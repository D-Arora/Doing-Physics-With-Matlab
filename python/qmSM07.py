# -*- coding: utf-8 -*-
"""
qmSM07.py    sep 2024

QUANTUM MECHANICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSM03.pdf


STATISTICAL MECHANICS: MB, BE, FD DISTRIBUTIONS: 
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
import random
import sys

tStart = time.time()


#%% VARIABLES
E = linspace(0.001,2,599)  # Energy grid: Emin to Emax  [eV]
Ef = 1.0                   # Fermi level [eV]
kB = 1.381e-23             # Boltzmann constant
e  = 1.602e-19             # electron charge

#%% MAXWELL - BOLTZMANN   Figure 1

def MB(T):
    fn = exp(-E*e/(kB*T))
    area = simps(fn,E)
    fn = fn/area
    return fn

plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('E',fontsize = 12)
ax.set_ylabel('f$_{MB}$',fontsize = 14)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([0,2])
#ax.set_xlim([0,nE])
T = 300; xP = E; yP = MB(T)
ax.plot(xP,yP,'b',lw = 2,label = T)
T = 3000; xP = E; yP = MB(T)
ax.plot(xP,yP,'r',lw = 2,label = T)
ax.legend()
fig1.tight_layout()


#%%  BOSE - EINSTEIN     Figure 2

def BE(T):
    fn = 1/(exp(E*e/(kB*T)) - 1)
    area = simps(fn,E)
    fn = fn/area
    return fn

plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('E',fontsize = 12)
ax.set_ylabel('f$_{BE}$',fontsize = 14)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([0,2])
#ax.set_xlim([0,nE])
T = 300; xP = E; yP = BE(T)
ax.plot(xP,yP,'b',lw = 2,label = T)
T = 3000; xP = E; yP = BE(T)
ax.plot(xP,yP,'r',lw = 2,label = T)
ax.legend()
fig2.tight_layout()


#%% FERMI - DIRAC    Figure 3

def FD(T):
    fn = 1/(exp((E-Ef)*e/(kB*T)) + 1)
    area = simps(fn,E)
    fn = fn/area
    return fn

plt.rcParams["figure.figsize"] = (5,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('E',fontsize = 12)
ax.set_ylabel('f$_{FD}$',fontsize = 14)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([0,2])
#ax.set_xlim([0,nE])
T = 300; xP = E; yP = FD(T)
ax.plot(xP,yP,'b',lw = 2,label = T)
T = 3000; xP = E; yP = FD(T)
ax.plot(xP,yP,'r',lw = 2,label = T)
ax.legend()
fig3.tight_layout()


#%% BM   BE   FD   distributions    Figure 4
# Prob E < 0.25  areaMB   areaBE    areaFD

target = 0.25
index1 = np.argmin(np.abs(E - target))
R = range(0,index1)
   
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('E',fontsize = 12)
ax.set_ylabel('f',fontsize = 14)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([0,2])
#ax.set_xlim([0,nE])
T = 3000; xP = E; yP = MB(T)
areaMB = simps(yP[R],E[R])
ax.plot(xP,yP,'b',lw = 2,label = 'M-B')

T = 3000; xP = E; yP = BE(T)
areaBE = simps(yP[R],E[R])
ax.plot(xP,yP,'r',lw = 2,label = 'B-E')

T = 3000; xP = E; yP = FD(T)
areaFD = simps(yP[R],E[R])
ax.plot(xP,yP,'k',lw = 2,label = 'F-D')

ax.legend()
fig4.tight_layout()

#%%  Console output
print('Probability E < 0.25 occupied')
print('   MB       BE       FD')
print('  %0.3f' % areaMB + '    %0.3f ' % areaBE + '   %0.3f'  % areaFD)


#%% SAVE FIGURES
#  fig1.savefig('a1.png')
#  fig2.savefig('a2.png')
#  fig3.savefig('a3.png')
#  fig4.savefig('a4.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)



