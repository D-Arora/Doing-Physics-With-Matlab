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


#%%
def fn(B):            # Calculation of total energy of a microstate
   # EB = E*B
    S = sum(E*B)
    return S

nP = 3               # Number of particles
num = 5               # Number of trials
nE = 4                # Number of energy levels
S = 0                   # Initial total energy of a microstate 
E = np.arange(0,nE,1)   # Energy levels
B = zeros([num,nE])     # Microstates
Emax = nE-1             # Total energy of a microstate

Es = 0; b = zeros(nE);b[0] = nP
for k in range(num):
    for n in range(0,nP):
        while Es != Emax:
           e = random.randint(0,nE-1)
           if e == 0:
             b[e] = b[e]
           if e > 0:
             b[e] = b[e]+1  
           Es = Es + e; 
          
           if Es > Emax:
              Es = Es - e; b[e] = b[e] - 1; 
        
        B[k,:] = b
        Es = 0; b = zeros(nE); b[0] = nP
        
# for k in range(num):
#     c = 0; b = zeros(nE) 
#     while c < nP:
#       e = random.randint(0,nE-1)
#       b[e] = b[e] + 1
#       S = fn(b)   
#       if S < Emax:
#          c = c+1
#       if S > Emax:
#          c = c-1 
#          b[e] = b[e] - 1
#       if S == Emax:
#          c = 100
#     B[k,:] = b
     
Bsum = np.sum(B,axis = 0)
Btotal = sum(Bsum)
prob = Bsum/Btotal

#%%
k = 3
K = 1
EB = linspace(0.1,50,599)
fBE = 1/(K*exp(k*EB)-1)


#%% Figure 1: Plot of Maxwell speed distribution
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('Energy E',fontsize = 12)
ax.set_ylabel('Probability prob(E)',fontsize = 12)
ax.set_ylim([0,0.2])
ax.set_xlim([0,nE])
for c in range(nE):
    xP = [c,c]
    yP = [0,prob[c]]
    ax.plot(xP,yP,'b',lw = 6)
ax.plot(EB,fBE,'r')
fig1.tight_layout()






#%% SAVE FIGURES
#  fig1.savefig('a1.png')
#  fig2.savefig('a2.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)



