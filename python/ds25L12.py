# -*- coding: utf-8 -*-
"""
ds25L12.py
Sep 25

DYNAMICAL SYSTEMS:
    Staddle node bifurcations: ghosts and bottlenecks


# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L12.pdf

"""

import numpy as np
from scipy.integrate import odeint, ode
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace, sqrt
from scipy.optimize import fsolve
from scipy import integrate

plt.close('all')

tStart = time.time()

#%%  SOLVE ODE x
def lorenz(t, state):    
    x = state
    dx = w - a*sin(x)
    return dx 

def XDOT(x):
    xDot = w - a*sin(x)
    return xDot

#%%  CELL 1  w > a
#  w > a   not fixed points, non uniform motion, period T
num = 599
w = 1
tMax = 500
a = 0.99
wa = w/a
x = linspace(0,2*pi,num)
xDot = XDOT(x)
A = linspace(0.1,0.9999,num)
T = np.zeros(num)
for c in range(num):
    fn = 1 / (w - A[c]*sin(x))
    T[c] = integrate.simpson(fn,x)

# Fig 1:a vs T
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('a')
ax.set_ylabel('T')
ax.set_title('$\omega$ = 1 ')
ax.grid()
ax.plot(A,T,'k',lw = 2)
fig1.tight_layout()
fig1.savefig('a1.png')

# Fig 2:x vs xDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.2)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x / $\pi$')
ax.set_ylabel('x$_{Dot}$ / $\pi$')
ax.set_title('$\omega$ = %0.3f' %w + '   a = %0.3f' %a 
             + '   $\omega$/a = %0.3f' %wa , fontsize = 12)
ax.grid()
ax.set_xlim([0,2])
ax.set_xticks(np.arange(0,2.2,0.25))
ax.plot(x/pi,xDot/pi,'k',lw = 2)
ax.plot([0,2],[0,0],'m',lw = 2)
fig2.tight_layout()
fig2.savefig('a2.png')

# Fig 3: t vs x
N = 9999
x0 = 0
t = linspace(0,tMax,N)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0]
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,2.2)
fig3, ax = plt.subplots(nrows=1, ncols=2)

plt.suptitle('$\omega$ = %0.3f' %w + '   a = %0.3f' %a 
             + '   $\omega$/a = %0.3f' %wa , fontsize = 12)
C = 0
ax[C].set_xlabel('t')
ax[C].set_ylabel('x / $\pi$')

ax[C].grid()

ax[C].plot(t,xS/pi,'k',lw = 2)
ax[C].plot([0,tMax],[0.5,0.5],'m',lw = 2)

C = 1
ax[C].set_xlabel('t')
ax[C].set_ylabel('x / $\pi$')
ax[C].grid()

ax[C].plot(t,(xS/pi) %(2),'k',lw = 2)
ax[C].plot([0,tMax],[0.5,0.5],'m',lw = 2)


fig3.tight_layout()
fig3.savefig('a3.png')

