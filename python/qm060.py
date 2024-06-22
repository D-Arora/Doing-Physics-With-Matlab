# -*- coding: utf-8 -*-
"""
qm050.py    June 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm060.pdf

QUANTUM MECHANICS     ROTATION SPECTRA: HCl molecule
  
"""

# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt
import scipy.special 
import scipy.special as sps
from scipy import *
import scipy.linalg as la
import math
import cmath
import os
import matplotlib.pyplot as plt
from pylab import *
from math import *
from numpy import linspace,matrix,array
from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve
from scipy.integrate import odeint, solve_ivp, simps
from scipy.sparse.linalg import eigsh, eigs #Solves the Eigenvalue problem
from scipy.sparse import diags #Allows us to construct our matrices
from matplotlib.animation import FuncAnimation, PillowWriter 
import time

tStart = time.time()

#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
cL = 2.99792458e8              # speed of light
amu = 1.67377e-27              # atomic mass unit  [kg]
se = e                         # Energy scaling factor   J <---> ev 
sx = 1e-9                      # Length scaling factor   m <---> nm

mH = 1.0078*amu                # mass H atom     [kg]
mCl = 34.9688*amu              # mass Cl atom   [kg]
R = 0.127e-9                   # Bond length

#%% SETUP 
mu = mH*mCl/(mH+mCl)           # reduced mass HCl molecule [kg]
I = mu*R**2                    # Moment of inertia
B = hbar**2/(2*I)


N = 25                        # Number of J values
J  = arange(1,N+1,1)
Js = arange(1,N+1,1)

Er = zeros(N)                 # Rotation energy levels   [J]
for q in range(N):
    Er[q] = q*(q+1)*B
    
ER = Er/se                   # Rotation energy levels    [eV]

# Emission spectra           J (intial state) -->  J-1 (final state
dE = (hbar**2/I)*Js/se
wL = 2*pi*I*cL/(hbar*Js)
k = 1/wL
f = cL/wL


#%%
EV = -array([4.39128216, 4.03387283, 3.67651618, 3.31921221, 2.96196096,
       2.60476241, 2.2476166 , 1.89052352, 1.53348319, 1.17649563,
       0.81956083, 0.46267883, 0.10584962])


#%%  FIG 1: ROTATION ENERGY LEVELS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, axes = plt.subplots(nrows=1, ncols=2)
C = 0   
xP = [0,1]
for q in range(12):
    yP = [ER[q],ER[q]]
    axes[C].plot(xP,yP, 'blue')
axes[C].set_xticks([])  
axes[C].set_ylabel('$E_R$  [ eV ]')       

C = 1
xP = J; yP = ER
axes[C].plot(xP,yP, 'bo',ms = 3)
axes[C].set_xlabel('J ')
axes[C].set_ylabel('$E_R$  [ eV ]')      
axes[C].set_xlim([0, 25])
axes[C].set_xticks(arange(0,27,5))
axes[C].xaxis.grid()
axes[C].yaxis.grid()
fig1.tight_layout()


# FIG 2:  ROTATIONAL ENERGY LEVEL SPACINGS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('dE  [ eV ]',color= 'black')
ax.set_xlabel('$J^*$',color = 'black')
ax.set_xlim([0,25])
ax.set_ylim([0,0.08])
xP = Js; yP = dE
ax.plot(xP,yP,'bo',ms = 3)
fig2.tight_layout()

# FIG 3:
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6.5,3)
fig3, axes = plt.subplots(nrows=1, ncols=2)
C = 0   
xP = Js; yP = wL*1e6
axes[C].plot(xP,yP, 'bo',ms = 3) 
axes[C].set_ylabel('$\lambda$  [ $\mu$m ]',color = 'b')  
axes[C].set_xlabel('$J^*$ ')
axes[C].set_xlim([0, 25]) 
axes[C].set_ylim([0, 500])  
axes[C].set_xticks(arange(0,26,5))
axes[C].xaxis.grid()
axes[C].yaxis.grid()  
axes[C].tick_params(axis='y', labelcolor = 'b')
axes[C].yaxis.grid(color = 'b')

ax2 = axes[C].twinx()
ax2.set_ylabel('f  [ THz ]',color= 'red')
# ax2.set_xlabel('t  [s]',color = 'red')
# ax2.set_xlim([0, 11])
# ax2.set_ylim([-70, 60])
# ax2.set_xticks(np.arange(0,11,2))
# ax2.set_yticks(np.arange(-80,61,20))
ax2.tick_params(axis='y', labelcolor = 'red')
ax2.yaxis.grid(color = 'r')
xP = Js; yP = f/1e12
ax2.plot(xP,yP,'ro',ms = 3)
fig3.tight_layout()

C = 1
xP = Js; yP = k/100
axes[C].plot(xP,yP, 'ko',ms = 3)
axes[C].set_xlabel('$J^*$ ')
axes[C].set_ylabel('k [ cm$^{-1}$ ]')      
axes[C].set_xlim([0, 25])
axes[C].set_xticks(arange(0,26,5))
axes[C].xaxis.grid()
axes[C].yaxis.grid()
fig3.tight_layout()


#%%  FIG 4:VIBRATION - ROTATION STATES
# Vibration
m = 3; n = 2

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig4, axes = plt.subplots(nrows=1, ncols=1)
xP = [0,1]
yP = [ -EV[m],-EV[m]]
axes.plot(xP,yP,'b',lw = 4)
yP =[ -EV[n],-EV[n]]
axes.plot(xP,yP,'b',lw = 4)
axes.set_xticks([]) 
axes.set_ylim([3.2,4]) 
axes.set_ylabel('$E_V$   $E_R$     [ eV ]')
axes.set_title('Vibration states:  m = %2.0f' % m + '   n = %2.0f' %n,  
               fontsize = 10, color = 'blue')
for q in range(15):
    yP = -EV[m] +[ ER[q],ER[q] ]
    axes.plot(xP,yP,'r',lw = 0.5)
    yP = -EV[n] +[ ER[q],ER[q] ]
    axes.plot(xP,yP,'m',lw = 0.5)
       

#%% TRANSITIONS BETWEEN STATES

# R branch
num = 10                      # number of rotation states
Ei = zeros(num); Ef = zeros(num)      
Ef = zeros(num)
Ei = EV[n] + ER[1:num+1]
for q in range(num):
     Ef[q] = EV[m] + ER[q]
Er = Ef - Ei
fr = Er*se / h

# P branch
Ei = zeros(num); Ef = zeros(num)      
Ef = zeros(num)
Ei = EV[n] + ER[0:num]
for q in range(num):
     Ef[q] = EV[m] + ER[q+1]
Ep = Ef - Ei
fp = Ep*se / h


# FIG 5:  VIBRATION-ROTATIONAL ABSORPTION SPECTRUM
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig5, ax = plt.subplots(1)

ax.xaxis.grid()
#ax.yaxis.grid()

ax.set_xlabel('$f_r$   $f_p$   [ x10$^{13}$  Hz] ',color = 'black')

#ax.set_xlim([8,9])
ax.set_ylim([0,1])
yP = [0,1]
for q in range(num):
    xP = [fr[q]/1e13,fr[q]/1e13]; 
    ax.plot(xP,yP,'r',lw = 3)
    xP = [fp[q]/1e13,fp[q]/1e13]; 
    ax.plot(xP,yP,'b',lw = 3)
ax.set_yticks([]) 
fig5.tight_layout()





#%%  SAVE FIGURES
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)

#%%
z = (fr[1]-fr[2])
print('z = %2.3e' % z)

z = hbar/(2*pi*I)
print('z = %2.3e' % z)

