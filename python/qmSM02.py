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

tStart = time.time()

#%%
# Number of energy values from E = 0 to Emax
num = 999
# Total energy of system of distinguishable non-interacting particles [eV]
Emin = 0; Emax = 0.4
# Boltzmann constant
kB = 1.381e-23
# electron charge
e = 1.602e-19
# Normalization constant (start value A = 1)
A = 1
# Total energy of system of particles [J]
E = e*linspace(Emin,Emax,num)
# Number of particles within unit volume V = 1
N = 1e8

# Function: Maxwell-Boltzmann distribution: fixed number of particles N
def fn(T):
    k = -1/(kB*T)
    f = exp(k*E)
    I = simps(f,E)
    A = N/I
    f = A*exp(k*E)
    return f

# Temperatures [K]
T = np.array([300,500,700,900])

# Maxwell-Boltzmann distribution for the 4 tempeartures
fMB = zeros([num,4])
for c in range(4):
    fMB[:,c] = fn(T[c])

# Number of particles in ranges E1 = 2e-20; E2 = 4e-20; E = emax = 6e-20
#    0 to E1   E1 to E2   E2 to E 
target_value = 2e-20
index1 = np.argmin(np.abs(E - target_value))
target_value = 4e-20
index2 = np.argmin(np.abs(E - target_value))
N1 = zeros(4); p =zeros(4)

for c in range(4):
    x1 = E[0:index1] 
    fn = fMB[0:index1,c]
    N1[c] = simps(fn,x1)
    p[c] = 100*N1[c]/N
print('Percentage particles 0 to E1')     
print('  T [K]          %0.0f' %T[0] + '    %0.0f' %T[1] +
           '    %0.0f' %T[2] + '   %0.f' %T[3] )
print('  p              %0.1f' %p[0] + '   %0.1f' %p[1] +
           '   %0.1f' %p[2] + '  %0.1f' %p[3] )   

for c in range(4):
    x1 = E[index1:index2] 
    fn = fMB[index1:index2,c]
    N1[c] = simps(fn,x1)
    p[c] = 100*N1[c]/N
print('Percentage particles E1 to E2')     
print('  T [K]          %0.0f' %T[0] + '    %0.0f' %T[1] +
           '    %0.0f' %T[2] + '   %0.f' %T[3] )
print('  p              %0.1f' %p[0] + '    %0.1f' %p[1] +
           '    %0.1f' %p[2] + '  %0.1f' %p[3] ) 

for c in range(4):
    x1 = E[index2:num] 
    fn = fMB[index2:num,c]
    N1[c] = simps(fn,x1)
    p[c] = 100*N1[c]/N
print('Percentage particles E2 to E')     
print('  T [K]          %0.0f' %T[0] + '    %0.0f' %T[1] +
           '    %0.0f' %T[2] + '   %0.f' %T[3] )
print('  p              %0.1f' %p[0] + '    %0.1f' %p[1] +
           '    %0.1f' %p[2] + '   %0.1f' %p[3] ) 
   
              
#%%  Figure 1: Plot of Maxwell-Boltzmann distribution
plt.rcParams["figure.figsize"] = (4,4)
fig1, ax = plt.subplots(nrows=1, ncols=1)
#fig1.subplots_adjust(top = 0.94, bottom = 0.15, left = 0.180,\
#                    right = 0.92, hspace = 0.36,wspace=0.40)

ax.set_xlabel('Energy E [ J ]',fontsize = 12)
ax.set_ylabel('particles/energy',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
xP = E
yP = fMB[:,0]
ax.plot(xP,yP,'b',lw = 2,label = '300')
yP = fMB[:,1]
ax.plot(xP,yP,'r',lw = 2,label = '500')
yP = fMB[:,2]
ax.plot(xP,yP,'m',lw = 2,label ='700')
yP = fMB[:,3]
ax.plot(xP,yP,'k',lw = 2,label = '900')
ax.legend()


fig1.tight_layout()

# fig1.savefig('a1.png')




#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)



