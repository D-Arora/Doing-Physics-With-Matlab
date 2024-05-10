# -*- coding: utf-8 -*-
"""
qm030.py            May 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

QUANTUM MECHANICS
   Scattering from a finite setp potential
   Numerical solutions
   
   https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm030.pdf

"""

# https://www.shivajicollege.ac.in/sPanel/uploads/econtent/247cb3d6399e1d4af61b59b5a12206c8.pdf


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



#%% FUNCTIONS
def lorenz(x, state):
      
    u, v = state
    P = E0
    
    if x > 0 and x < aB:
       P = E0 - U0
   
    du = v 
    dv = C*P*u
    return [du, dv]



#%%  INPUTS
# X domain  [nm]
xMin = -1.0; xMax = 1.0 
# beam particle energy  [eV]   
E0 = 5  
# height of potential barrier [eV]         
U0 = -50  
# Barrier width [nm]         
aB = 1.2*0.2676280855685051    
# aB = 0.2
# Grid points
N = 18129                           

#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [ J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
me = 9.10938356e-31            # Electron mass [kg]
C = -2*me*e/hbar**2            # SE constant 
           
#%% COMPUTATIONS
xMin = xMin*1e-9; xMax = xMax*1e-9    # [nm --> m]
x = np.linspace(xMin,xMax,N)
aB = aB*1e-9                          # [nm --> m]
E = e*E0                          # beam particle energy  [eV -->J]                                 
U = e*U0                          # height of potential step [eV -->J]

w = E/hbar                        # angular frequency [rad/s]
P = 2*pi/w                        # inc wave: oscillation period [s]
k1 = sqrt(2*me*E)/hbar            # propagation constant: regions 1 & 3 [1/m]
k2 = cmath.sqrt(2*me*(E-U))/hbar  # propagation constant: region 2     [1/m]
L1 = 2*pi/k1                      # wavelength: regions 1 & 3    [m]
L2 = 2*pi/k2                      # wavelength: region 2         [m]    


#%% SOLVE SCHRODINGER EQUATION
xSpan = np.flip(x)
# Real part
u1 = 0; u2 = k1
u0 = [u1,u2]
sol = odeint(lorenz, u0, xSpan,  tfirst=True)
psiR = sol[:,0]
# Imaginary part
u1 = 1; u2 = 0
u0 = [u1,u2]
sol = odeint(lorenz, u0, xSpan, tfirst=True)
psiI = sol[:,0]

psiR = np.flip(psiR)  
psiI = np.flip(psiI) 

probD = psiR**2 + psiI**2               # Probability density


# Probability per unit length
fn = probD[x<0]
z = x[x<0]
prob1 = simps(fn,z)
fn = probD[x>aB]
z = x[x>aB]
prob2 = simps(fn,z)
T21 = 100*(prob2/prob1)*(-xMin/(xMax-aB))   # % Transmission probability 

#%% CONSOLE DISPLAY
print('E0 = %2.0f  eV' % E0)
print('U0 = %2.0f  eV' % U0)
s = aB*1e9; print('a  = %2.3f  nm' % s)
s = L1*1e9; print('wavelength: regions 1 & 3  wL = %2.3f  nm' %s)
if E0 > U0:
      s = real(L2)*1e9; print('wavelength: region 2  wL2 = %2.3f nm' % s)
      s = aB/real(L2);  print('a/wL2  = %2.3f ' % s)
print('Transmission probability percentage T = %2.0f' % T21)

#%% GRAPHICS
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6,4)

fig1, axes = plt.subplots(nrows=2, ncols=1)
xP = x*1e9

R = 0
axes[R].set_ylabel('$\psi$',color = 'black',fontsize = 16)
axes[R].xaxis.grid()
axes[R].yaxis.grid()
axes[R].set_xlim(xMin*1e9,xMax*1e9)
axes[R].set_xticks(np.arange(-1.0,1.1,0.2))
yP = psiR
axes[R].plot(xP, yP,'b', lw = 2)
yP = psiI
axes[R].plot(xP, yP,'r', lw = 2)
xx = aB*1e9; yy = max(psiI)
axes[R].plot([xx,xx],[-yy,yy],'k',lw = 1)
axes[R].plot([0,0],[-yy,yy],'k',lw = 1)
R = 1
axes[R].set_xlabel('x  [ nm ] ',color = 'black',fontsize = 12)
axes[R].set_ylabel('|$\psi|^2$',color = 'black',fontsize = 14)
axes[R].xaxis.grid()
axes[R].yaxis.grid()
axes[R].set_xlim(xMin*1e9,xMax*1e9)
yP = probD
axes[R].plot(xP, yP,'b', lw = 2)
axes[R].fill_between(xP, probD,color = [1,0,1],alpha=0.2)
axes[R].plot([0,0],[0,max(probD)],'k',lw = 1)
xP = aB*1e9
axes[R].plot([xP,xP],[0,max(probD)],'k',lw = 1)
fig1.tight_layout()

fig1.savefig('a1.png')


#%%
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)

fig1, axes = plt.subplots(nrows=1, ncols=1)
axes.set_ylabel('T, R',color = 'black',fontsize = 14)
axes.set_xlabel('E / U',color = 'black',fontsize = 14)
axes.xaxis.grid()
axes.yaxis.grid()

num = 599
U1 = 10
a = 10e-10

EA = 0; EB = U1-0.01
E1 = linspace(EA,EB,num)

T1 = sqrt(2*me*U1*e*a**2/hbar**2)
T2 = sqrt(1-E1/U1)
ka = T1*T2
T3 = np.sinh(ka)**2
T4 = 4*(E1/U1)*(1-E1/U1)
T5 = 1+T3/(T4+1e-16)
T = 1/(T5+1e-16)
R = 1 - T
axes.plot(E1/U1,T,'b',lw = 2, label ='T')
axes.plot(E1/U1,R,'r',lw = 1, label = 'R')
axes.legend()

#%%
EA = U1+0.1; EB = 5*EA
E1 = linspace(EA,EB,num)
T1 = sqrt(2*me*U1*e*a**2/hbar**2)
T2 = sqrt(E1/U1-1)
ka = T1*T2
T3 = sin(ka)**2
T4 = 4*(E1/U1)*(E1/U1-1)
T5 = 1+T3/(T4+1e-16)
T = 1/T5
R = 1 - T
aS = a*1e9
axes.plot(E1/U1,T,'b',lw = 2)
axes.plot(E1/U1,R,'r',lw = 1)
axes.set_title('a = %2.2f nm' % aS + '     $U_0$ = %2.0f eV' % U1,
               fontsize = 12, color = 'black')
fig1.tight_layout()

fig1.savefig('a2.png')
    

    
#%%
# z = ones(7)
# z[5]=5
# zz = zeros(9)
# zz[1:8]=z
# print(z)
# print(zz)


