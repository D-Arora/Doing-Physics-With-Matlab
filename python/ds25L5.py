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

plt.close('all')

tStart = time.time()


#%%  SOLVE ODE x
def lorenz(t, state):    
    x = state
    dx = x-x**3
    return dx


#%% SETUP
u0 = 2
tMax = 5; N = 9999

#%% SOLVE ODE
t = linspace(0,tMax,N)
sol = odeint(lorenz, u0, t, tfirst=True)
x = sol[:,0] 

# fixed points
xss = zeros(3)
xss[0] = 1; xss[1] = -1; xss[2] = 0

# xDot
X = linspace(-2, 2,999)
Xdot = X - X**3

#xxx
#%% FIGURE 1: t vs x
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('t')
ax.set_ylabel('x')
ax.grid()

ax.plot(t,x,'k',lw = 2)

xP = [0,tMax]; yP = [xss[0], xss[0]]; ax.plot(xP,yP,'b',lw = 1)
yP = [xss[1], xss[1]]; ax.plot(xP,yP,'b',lw = 1)
yP = [xss[2], xss[2]]; ax.plot(xP,yP,'r',lw = 1)
fig1.tight_layout()

#xxx
#%% FIGURE 2: x vs xDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('x')
ax.set_ylabel('xDot')
ax.grid()
ax.plot(X,Xdot,'k',lw = 2)

ax.plot(xss[0],0,'bo',ms = 8)
ax.plot(xss[1],0,'bo',ms = 8)
ax.plot(xss[2],0,'ro',ms = 8)
ax.arrow(-2, 0.9, 0.3, 0,  lw = 3,head_width=0.4, head_length=0.2, fc='k', ec='k')
ax.arrow(0.2, 0.9,0.3, 0,  lw = 3,head_width=0.4, head_length=0.2, fc='k', ec='k')
ax.arrow(-0.2, 0.9, -0.3, 0,  lw = 3,head_width=0.4, head_length=0.2, fc='k', ec='k')
ax.arrow(2, 0.9, -0.3, 0,   lw = 3,head_width=0.4, head_length=0.2, fc='k', ec='k')
fig2.tight_layout()

#%% POTENTIAL
V = -(X**2/2 - X**4/4)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('V')
ax.grid()

ax.plot(X,V,'k',lw = 2)
yP = -(xss[0]**2/2 - xss[0]**4/4); ax.plot(xss[0],yP,'bo',ms = 8)
yP = -(xss[1]**2/2 - xss[1]**4/4); ax.plot(xss[1],yP,'bo',ms = 8)
yP = -(xss[2]**2/2 - xss[2]**4/4); ax.plot(xss[2],yP,'ro',ms = 8)
fig3.tight_layout()


#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
