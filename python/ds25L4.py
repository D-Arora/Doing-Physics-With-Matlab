# -*- coding: utf-8 -*-
"""
ds25L4.py
Aug 25

DYNAMICAL SYSTEMS:
    


# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L4.pdf

"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, tan, linspace, sqrt, zeros
from scipy.optimize import fsolve

plt.close('all')

tStart = time.time()



#%% SETUP
tMax = 20; N = 999
t = linspace(0,tMax,N)
xS = (2*t/3)**(3/2) 

x = linspace(0.01,400,N)
xDot = x**(1/3)

dfdx = 1/(3*x**(2/3))
    
#%% FIGURE 1: t vs x
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.grid()
ax.plot(t,xS,'k',lw = 2)
fig1.tight_layout()

#%% FIGURE 2: x vs xDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('x$_{dot}$')
ax.grid()
ax.plot(x,xDot,'k',lw = 2)
fig2.tight_layout()

#%% FIGURE 3:
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('df / dx')
ax.grid()
ax.plot(x,dfdx,'k',lw = 2)
fig3.tight_layout()


#%%
from sympy import symbols, integrate

x = symbols('x')
eq = 1/(1+x**2)
integral = integrate(eq, x)

t = linspace(-0.98*pi/2,0.98*pi/2,N)
xt = tan(t)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t / $\pi$')
ax.set_ylabel('x ')
ax.set_xticks(np.arange(-0.5,0.52,0.25))
ax.grid()
ax.plot(t/pi,xt,'k',lw = 2)
fig4.tight_layout()

#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
