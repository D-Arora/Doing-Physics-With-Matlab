# -*- coding: utf-8 -*-
"""
qmSM01.py    Aug 2024

QUANTUM MECHANICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSM01.pdf


STATISTICAL MECHANICS: MAXWELL - BOLTZMANN DISTRIBUTION
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
# Total number of particles
N = 6
# Number of energy levels: energy levels 0 1 2 3 4 5 6 7 8
nE = 9
# Number of macroStates
nMacro = 20

# Populations of macrostates: number of particles for each energy level
macroS = zeros([nMacro,nE])
macroS[0,:]  = np.array([5,0,0,0,0,0,0,0,1])
macroS[1,:]  = np.array([4,1,0,0,0,0,0,1,0])
macroS[2,:]  = np.array([4,0,1,0,0,0,1,0,0])
macroS[3,:]  = np.array([3,2,0,0,0,0,1,0,0])
macroS[4,:]  = np.array([4,0,0,1,0,1,0,0,0])
macroS[5,:]  = np.array([3,1,1,0,0,1,0,0,0])
macroS[6,:]  = np.array([2,3,0,0,0,1,0,0,0])
macroS[7,:]  = np.array([4,0,0,0,0,2,0,0,0])
macroS[8,:]  = np.array([3,1,0,1,1,0,0,0,0])
macroS[9,:]  = np.array([3,0,2,0,1,0,0,0,0])
macroS[10,:] = np.array([2,2,1,0,1,0,0,0,0])
macroS[11,:] = np.array([1,4,1,0,1,0,0,0,0])
macroS[12,:] = np.array([3,0,1,2,0,0,0,0,0])
macroS[13,:] = np.array([2,2,0,2,0,0,0,0,0])
macroS[14,:] = np.array([2,1,2,1,0,0,0,0,0])
macroS[15,:] = np.array([1,3,1,1,0,0,0,0,0])
macroS[16,:] = np.array([0,5,0,1,0,0,0,0,0])
macroS[17,:] = np.array([2,0,4,0,0,0,0,0,0])
macroS[18,:] = np.array([1,2,3,0,0,0,0,0,0])
macroS[19,:] = np.array([0,4,2,0,0,0,0,0,0])

# Number of microstates for each macrostate
microS = zeros(nMacro)
Nf = factorial(N)
nf = zeros(9)
for c1 in range(nMacro):
    q = 1
    for c2 in range(nE):
       q = q*factorial(macroS[c1,c2]) 
    microS[c1] = Nf / q          
# total number of microstates
microN = sum(microS)

# Probability of observing of a given macrostate
p = zeros(nMacro)
nAvg = zeros(nE)
for c in range(nMacro):
    p[c] = microS[c]/microN

# Average number of particles in a given macrostate
for c in range(nE):
    nAvg[c] = sum(macroS[:,c]*p)

# Probability of finding a particle with a given energy
probE = nAvg/6

# Energy levels
EL = np.arange(0,9,1)

E = linspace(0,8,299)
A = probE[0]
k1 = -np.log(probE[8]/A)/8
k2 = -np.log(probE[2]/A)/2
k3 = -np.log(probE[1]/A)/1
k = (k1+k2+k3)/3
y = A*exp(-k*E)

# Console output
print('Number of microstates for each energy level (macrostate)')
print(microS)
print('Total number of microstates')
print(microN)
print('  ')
print('Probability of finding a particle with a given energy') 
print(probE)


#%%  Figure 1: plot of Maxwell-Boltzmann distribution
plt.rcParams["figure.figsize"] = (4,4)
fig1, ax = plt.subplots(nrows=1, ncols=1)
#fig1.subplots_adjust(top = 0.94, bottom = 0.15, left = 0.180,\
#                    right = 0.92, hspace = 0.36,wspace=0.40)

ax.set_xlabel('energy',fontsize = 12)
ax.set_ylabel('Prob. partilce with energy E',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
xP = EL; yP = probE
ax.plot(xP,yP,'bo')
xP = E; yP = y
ax.plot(xP,yP,'r')
fig1.tight_layout()


# fig1.savefig('a1.png')





#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
