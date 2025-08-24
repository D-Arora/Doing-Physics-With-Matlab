# -*- coding: utf-8 -*-
"""
ds25L3.py
Aug 25

DYNAMICAL SYSTEMS:
    Fixed Points and Stability 


# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L3.pdf

"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace, sqrt, zeros
from scipy.optimize import fsolve

plt.close('all')

tStart = time.time()


#%%  SOLVE ODE x
def lorenz(t, state):    
    x = state
    dx = x-cos(x)  
    return dx


#%% SETUP
u0 = 0.73
tMax = 10; N = 9999

#%% SOLVE ODE
t = linspace(0,tMax,N)
sol = odeint(lorenz, u0, t, tfirst=True)
x = sol[:,0] 

# fixed points
def equations(variables):
    Z = variables  # Unpack the variables
    eq = Z - cos(Z)    
    return eq
IC = [1.0] # Initial guess for x and y
xss = fsolve(equations, IC)
print(xss)

# xDot
X = linspace(-pi, pi,999)
Xdot = X-cos(X)

#xxx
#%% FIGURE 1: t vs x
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('t')
ax.set_ylabel('x')
ax.grid()

ax.plot(t,x,'k',lw = 2)

# xP = [0,tMax]; yP = [xss[0], xss[0]]; ax.plot(xP,yP,'r',lw = 1)
# yP = [xss[1], xss[1]]; ax.plot(xP,yP,'b',lw = 1)

fig1.tight_layout()

#xxx
#%% FIGURE 2: x vs xDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('x  [rad]')
ax.set_ylabel('xDot  [rad/s]')
ax.set_title('x$_{ss}$ = %0.4f rad' % xss[0])
ax.grid()

ax.plot(X,X,'m',lw = 2)
ax.plot(X,cos(X),'g',lw = 2)
ax.plot(X,Xdot,'k',lw = 2)

xP = xss; yP = xss-cos(xss)
ax.plot(xP,yP,'ro',ms = 6)
ax.arrow(-0.5, 0, -0.3, 0,  lw = 3,head_width=0.4, head_length=0.2, fc='k', ec='k')
ax.arrow(2, 0, 0.3, 0,  lw = 3,head_width=0.4, head_length=0.2, fc='k', ec='k')


fig2.tight_layout()


#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
