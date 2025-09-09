# -*- coding: utf-8 -*-
"""
ds25L12N.py
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
    dx = r + x**2
    return dx 

def XDOT(x):
    xDot = r + x**2
    return xDot

#%%    
num = 599
r = 0.001
tMax = 49.6

x2 = 10
x = linspace(0,x2,num)
xDot = XDOT(x)

xT = linspace(0,10,num)
rT = linspace(0.001,10,num)
T = np.zeros(num)
for c in range(num):
    fn = 1 / (rT[c] + xT**2)
    T[c] = integrate.simpson(fn,x)

# Fig 1:a vs T
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r')
ax.set_ylabel('T')
ax.grid()
ax.plot(rT,T,'k',lw = 2)
fig1.tight_layout()
fig1.savefig('a1.png')

# Fig 2:x vs xDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.2)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('x$_{Dot}$')
ax.set_title('r = %0.3f' %r , fontsize = 12)
ax.grid()
# ax.set_xlim([0,2])
# ax.set_xticks(np.arange(0,2.2,0.25))
ax.plot(x,xDot,'k',lw = 2)
ax.plot([0,x2],[0,0],'m',lw = 2)
fig2.tight_layout()
fig2.savefig('a2.png')

#Fig 3: t vs x
N = 9999
x0 = 0
t = linspace(0,tMax,N)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0]
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,2.2)
fig3, ax = plt.subplots(nrows=1, ncols=1)

plt.suptitle('r = %0.3f' % r, fontsize = 12)

ax.set_xlabel('t')
ax.set_ylabel('x')
ax.grid()
ax.plot(t,xS/pi,'k',lw = 2)
#ax[C].plot([0,tMax],[0.5,0.5],'m',lw = 2)

fig3.tight_layout()
fig3.savefig('a3.png')

