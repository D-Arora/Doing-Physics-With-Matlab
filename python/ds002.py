# -*- coding: utf-8 -*-
"""
Created on Sun Aug 10 10:16:11 2025

@author: Owner
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace
from scipy.optimize import fsolve
plt.close('all')

tStart = time.time()


#%%  SOLVE ODE x
def lorenz(t, state):    
    x = state
    dx = w - a*sin(x)
    return dx 
#%% SETUP
x0 = 0
w = 2*pi
a = 0
tMax = 2

#%% SOLVE ODE
t = linspace(0,tMax,9999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
xS = xS % (2*pi)

# FIGURE 1: time evolution of xS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x / $\pi$')
ax.grid()
ax.plot(t,xS/pi,'b',lw = 2)
fig1.tight_layout()


#%% ODE
num = 999
X = linspace(0,2*pi,num)
Xdot = w - a*sin(X)

k = np.zeros(2); p = 0
for c in range(num-2):
    q = Xdot[c]*Xdot[c+1]
    if q <= 0:
       k[p] = c
       p = int(p+1)     
k = k.astype(int)
vZeros = X[k]    
print(vZeros/pi)
k
# FIGURE 1: X vs Xdot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x/$\pi$')
ax.set_ylabel('xDot/$\pi$')
ax.grid()
ax.plot(X/pi,Xdot/pi,'b',lw = 2)
fig2.tight_layout()





