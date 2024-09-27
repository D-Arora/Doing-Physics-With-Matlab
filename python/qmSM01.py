# -*- coding: utf-8 -*-
"""
qmSM01.py    sep 2024

QUANTUM MECHANICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSM03.pdf


STATISTICAL MECHANICS:
    MB   BE   FD DISTRIBUTIONS
    6 particles, 9 energy levels, Emax = 8
    Computation of probability that an energy state is occupied
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

#%%  MACROSTATES
# N       Total number of particles
# nE      Number of energy states: energy states (levels) 0 1 2 3 4 5 6 7 8
# nMacro  Number of macroStates
#macroS   # Populations of macrostates: number of particles for each energy state
N = 6
nE = 9
nMacro = 20

macroS = zeros([nMacro,nE])
macroS[0,:]  = np.array([5,0,0,0,0,0,0,0,1])
macroS[1,:]  = np.array([4,1,0,0,0,0,0,1,0])
macroS[2,:]  = np.array([4,0,1,0,0,0,1,0,0])
macroS[3,:]  = np.array([3,2,0,0,0,0,1,0,0])
macroS[4,:]  = np.array([4,0,0,1,0,1,0,0,0])
macroS[5,:]  = np.array([3,1,1,0,0,1,0,0,0])
macroS[6,:]  = np.array([2,3,0,0,0,1,0,0,0])
macroS[7,:]  = np.array([4,0,0,0,2,0,0,0,0])
macroS[8,:]  = np.array([3,1,0,1,1,0,0,0,0])
macroS[9,:]  = np.array([3,0,2,0,1,0,0,0,0])
macroS[10,:] = np.array([2,2,1,0,1,0,0,0,0])
macroS[11,:] = np.array([1,4,0,0,1,0,0,0,0])
macroS[12,:] = np.array([3,0,1,2,0,0,0,0,0])
macroS[13,:] = np.array([2,2,0,2,0,0,0,0,0])
macroS[14,:] = np.array([2,1,2,1,0,0,0,0,0])
macroS[15,:] = np.array([1,3,1,1,0,0,0,0,0])
macroS[16,:] = np.array([0,5,0,1,0,0,0,0,0])
macroS[17,:] = np.array([2,0,4,0,0,0,0,0,0])
macroS[18,:] = np.array([1,2,3,0,0,0,0,0,0])
macroS[19,:] = np.array([0,4,2,0,0,0,0,0,0])


#%% MAXWELL-BOLTZMANN STATISTICS
# microS     Number of microstates for each macrostate
# Nf         factorial(N)
# c          macrostates 0 to 19
# c2         energy states 0 to 8 
# q          number of micostates for each macrostate
# Nmicro     total number of microstates 
# probS      probability of obsering a macrostate                
# nAvg       average number of particles in a given macrostate for each energy state
# probMB     probability of finding a particle with a given energy

microS = zeros(nMacro)
Nf = factorial(N)

# Compute: number of microstates for each macrostate
for c in range(nMacro):
    q = 1
    for c2 in range(nE):
       q = q*factorial(int(macroS[c,c2])) 
    microS[c] = Nf / q          
Nmicro = sum(microS)

# Compute: Probability of observing of a given macrostate
probS = zeros(nMacro)
for c in range(nMacro):
    probS[c] = microS[c]/Nmicro

# Compute:Avg number of particles in a given macrostate for each energy state
nMB = zeros(nE)
for c in range(nE):
    nMB[c] = sum(macroS[:,c]*probS)

# Compute: probability of finding a particle with a given energy
probMB = nMB/N


#%% BOSE-EINSTEIN STATISTICS: probability distribution
# macroS    macrostates
# nBE       number of particles in each energy state
# num       total number of particles in each energy states
# probBE    prob of an energy state being occupied

nBE = np.sum(macroS,axis = 0)
num = np.sum(macroS, axis = None)
probBE = nBE/num


#%% FERMI-DIRAC STATISTICS: probability distribution
# macroFD     macrostates 
# nFD         number of particles in each energy state
# num         total number of particles in each energy states
# nFD         
macroFD = zeros([3,nE])
macroFD[0,:] = macroS[10,:]
macroFD[1,:] = macroS[13,:]
macroFD[2,:] = macroS[14,:]
nFD = np.sum(macroFD,axis = 0)
num = 18
probFD = nFD/num


#%% Compare probability distributions: relative probabilities
probMax = max(probBE)
probMB_R = probMB/probMax
probBE_R = probBE/probMax
probFD_R = probFD/probMax

EL = range(0,9,1)


#%%  Figure 1: plot of Maxwell-Boltzmann distribution
plt.rcParams["figure.figsize"] = (6,4)
fig1, ax = plt.subplots(nrows=1, ncols=1)
#fig1.subplots_adjust(top = 0.94, bottom = 0.15, left = 0.180,\
#                    right = 0.92, hspace = 0.36,wspace=0.40)

ax.set_xlabel('energy',fontsize = 12)
ax.set_ylabel('relative probabilities',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
xP = EL; yP = probMB_R
ax.plot(xP,yP,'bo')
xP = EL; yP = probBE_R
ax.plot(xP,yP,'ro')
xP = EL; yP = probFD_R
ax.plot(xP,yP,'ko')


# Fit cubicSoline to data
from scipy.interpolate import CubicSpline

cs = CubicSpline(xP, probMB_R)
x_range = np.arange(0, 9, 1)
plt.plot(x_range, cs(x_range),'b',label = 'MB')

cs = CubicSpline(xP, probBE_R)
x_range = np.arange(0, 9, 1)
plt.plot(x_range, cs(x_range),'r',label = 'BE')

cs = CubicSpline(xP, probFD_R)
x_range = np.arange(0, 9, 1)
ax.plot(x_range, cs(x_range),'k', label = 'FD')

ax.legend()

fig1.tight_layout()


#%% Console output
# print('Number of microstates for each energy level (macrostate)')
print('state     probMB_R     probBE_R     probFD_R') 
for c in range(9):
    print('  %0.0f' % c + '         %0.3f' % probMB_R[c] + 
          '        %0.3f' % probBE_R[c] + '       %0.3f'  % probFD_R[c])


#%% SAVE FIGURES
#  fig1.savefig('a1.png')
#  fig2.savefig('a2.png')
#  fig3.savefig('a3.png')
#  fig4.savefig('a4.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


