# -*- coding: utf-8 -*-
"""
qmSM08.py    sep2024

QUANTUM MECHANICS  /  STATISTICAL MECHANICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSM08.pdf


Specific heat of solids
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


#%%  SPECIFIC HEAT OF SOLIDS
R  = 8.3141         # Universal gas constant  [J/(mol.K)]
kB = 1.30805e-23   # Boltzmaan constant  [J/K]


# Debye( (Einstein) temperature TD:  x = T / TD   z = TD / T
N = 9999; xMin = 1e-3; xMax = 5
x = linspace(xMin,xMax,N)
z = 1/x

# Classical Model
CL = 234*R*x**3

# Einstein model
Cv = 3*R*(z)**2*exp(z) / (exp(z) - 1)**2

# Debye model
c1 = 36*R; c2 = 9*R; CV = zeros(N)
for c in range(N):
    uMax = z[c]
    u = linspace(xMin,uMax,999)
    fn =  u**3/(exp(u)-1)
    I = simps(fn,u)
    CV[c] = c1*x[c]**3 * I - c2*z[c]/(exp(z[c])-1)

# Metals
#theta = 165
#T = theta*x


#%% GRAPHICS
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('T/T$_D$',fontsize = 12)
ax.set_ylabel('C$_v$ / 3R',fontsize = 12)
ax.set_ylim([0,1.2])
ax.grid()
xP = x; yP = ones(N)
ax.plot(xP,yP,'m',lw = 1,label = 'C')
yP = Cv/(3*R)
ax.plot(xP,yP,'b',lw = 2,label = 'E')
yP = CV/(3*R)
ax.plot(xP,yP,'r',lw = 2, label ='D')
ax.legend()
fig1.tight_layout()


plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('T  [K]',fontsize = 12)
ax.set_ylabel('C$_v$ / 3R',fontsize = 12)
ax.set_ylim([0,1.2])
ax.set_xlim([0,300])
ax.grid()
theta = 428  # aluminum
T = theta*x
xP = T; yP = CV/(3*R)
ax.plot(xP,yP,'b',lw = 2, label ='Al')
theta = 165  # gold
T = theta*x
xP = T; yP = CV/(3*R)
ax.plot(xP,yP,'r',lw = 2, label ='Au')
theta = 105    # laed
T = theta*x
xP = T; yP = CV/(3*R)
ax.plot(xP,yP,'m',lw = 2, label ='Pb')

ax.legend()
fig2.tight_layout()

#%%
plt.rcParams["figure.figsize"] = (5,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('T/T$_D$',fontsize = 12)
ax.set_ylabel('C$_v$ / 3R',fontsize = 12)
ax.set_ylim([0,0.2])
ax.set_xlim([0,0.1])
ax.grid()
xP = x; yP = CV/(3*R)
ax.plot(xP,yP,'b',lw = 2, label ='Debye')

yP = 75*x**3
ax.plot(xP,yP,'r',lw = 2, label ='$T^3$')

ax.legend()
fig3.tight_layout()




#%% SAVE FIGURES
#  fig1.savefig('a1.png')
#  fig2.savefig('a2.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)



