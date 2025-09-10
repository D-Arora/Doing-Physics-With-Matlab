# -*- coding: utf-8 -*-
"""
ds25L13.py
Ian Cooper
Sep 25

DYNAMICAL SYSTEMS:
    Modelling Firefly Entrainment - Dynamical Systems

https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation

https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L11.pdf

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

#%% INPUTS
N = 9999
num = 999
wF = 1
wS = 0.5
A = 1
tMax = 10
x0 = 0.5*pi

#%% Fig. 1: x vs xDot
x = linspace(-pi,pi,num)
xDot = wS - wF - A*sin(x)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x / $\pi$')
ax.set_ylabel('x$_{Dot}$ / $\pi$')
ax.set_title('w$_S$ = %0.3f' %wS) 
ax.grid()
ax.plot(x/pi,xDot/pi,'k',lw = 2)
ax.plot([-1,1],[0,0],'m',lw = 2)
fig1.tight_layout()
fig1.savefig('a1.png')

#%% Fig. 2: solve ODe  t vs x

t = linspace(0,tMax,N)
def lorenz(t, state):    
    x = state
    dx = wS - wF - A*sin(x)
    return dx 

sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0]

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_ylabel('x/$\pi$')
ax.set_xlabel('t')
ax.set_title('w$_S$ = %0.3f' %wS) 
ax.grid()
ax.plot(t,xS/pi,'k',lw = 2)
ax.plot([0,tMax],[0,0],'m',lw = 2)
fig2.tight_layout()
fig2.savefig('a2.png')


