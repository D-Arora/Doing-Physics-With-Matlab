# -*- coding: utf-8 -*-
"""
qmSM06.py    Aug 2024

QUANTUM MECHANICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSM03.pdf


STATISTICAL MECHANICS: BOSE-EINSTEIN DISTRIBUTION: photons or phonons
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
T = 5000                   # Temperature [K]
E = linspace(0.001,3,599)  # Energy grid: Emin to Emax  [eV]
Ef = 1.0                   # Fermi level [eV]
kB = 1.381e-23             # Boltzmann constant
e  = 1.602e-19             # electron charge

#%% DISTRIBUTIONS
# Bose-Einstein
fBE = 1/(exp(E*e/(kB*T)) - 1)
areaBE = simps(fBE,E)
fBE = fBE/areaBE

# Fermi-Dirac
fFD = 1/(exp((E-Ef)*e/(kB*T)) + 1)
areaFD = simps(fFD,E)
fFD = fFD/areaFD

# Maxwell-Boltzmann
fMB = exp(-E*e/(kB*T))
areaMB = simps(fMB,E)
fMB = fMB/areaMB


#%% Figure 1: Plot of Bose-Einstein distribution
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('E',fontsize = 12)
ax.set_ylabel('f$_{BE}$',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([0,2])
#ax.set_xlim([0,nE])
xP = E; yP = fBE
ax.plot(xP,yP,'b',lw = 2)
fig1.tight_layout()


#%% Figure 2: Plot of Fermi-dirac distribution
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('E',fontsize = 12)
ax.set_ylabel('f$_{FD}$',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([0,1.00])
xP = E; yP = fFD
ax.plot(xP,yP,'b',lw = 2)
xP = [Ef,Ef]; yP = [0,0.5]
ax.plot(xP,yP,'r',lw = 1)
xP = [0,Ef]; yP = [0.5,0.5]
ax.plot(xP,yP,'r',lw = 1)
fig2.tight_layout()

#%% Figure 3: Plot of distributions
plt.rcParams["figure.figsize"] = (5,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('E',fontsize = 12)
ax.set_ylabel('prob densities, f',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([0,2])
#ax.set_xlim([0,nE])
xP = E; yP = fBE
ax.plot(xP,yP,'b',lw = 2,label='f$_{BE}$')
yP = fFD
ax.plot(xP,yP,'r',lw = 2,label='f$_{FD}$')
yP = fMB
ax.plot(xP,yP,'m',lw = 2,label='f$_{MD}$')
ax.legend()
fig3.tight_layout()


#%% FERMI-DIRAC simple model 6 balls and 9 energy states

s = np.arange(0,9,1)
probFD = zeros(9)
probFD[0:5] = np.array([6,5,3,3,1])/18

# Figure 4
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('E',fontsize = 12)
ax.set_ylabel('prob densities, f$_{FD}$',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([-0.02,0.4])
#ax.set_xlim([0,nE])
xP = s; yP = probFD
ax.plot(xP,yP,'ob',ms = 6,)

fig4.tight_layout()


#%% SAVE FIGURES
#  fig1.savefig('a1.png')
#  fig2.savefig('a2.png')
#  fig3.savefig('a3.png')
#  fig4.savefig('a4.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)



