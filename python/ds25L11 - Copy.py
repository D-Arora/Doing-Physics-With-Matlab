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
from scipy import integrate

plt.close('all')

tStart = time.time()


#%%  SOLVE ODE x
def lorenz(t, state):    
    x = state
    dx = w - a*sin(x)
    return dx 
#%% SETUP
x0 = 1*pi
tMax = 20; N = 9999
w = 1
a = 0.9999

col = [0,0,1]

#%% SOLVE ODE
t = linspace(0,tMax,N)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
#xS = xS % (2*pi)

#xxx
# FIGURE 1: time evolution of xS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x / $\pi$')
ax.grid()

# ax.set_ylim([0,0.52])
# ax.set_yticks(np.arange(0,0.52,0.1))

ax.plot(t,xS/pi,color = col,lw = 2)

fig1.tight_layout()


#xxx

#%% Find fixed points
num = 9999
X = linspace(0,2*pi,num)
Xdot = w - a*sin(X)
#Xdot = X
k = np.zeros(2); p = 0
for c in range(num-2):
    q = Xdot[c]*Xdot[c+1]
    if q <= 0:
       k[p] = c
       p = int(p+1)     
k = k.astype(int)
Z = X[k]    

#%%

      
# FIGURE 2: X vs Xdot
plt.rcParams['font.size'] = 14
plt.rcParams["figure.figsize"] = (5,2.2)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x / $\pi$')
ax.set_ylabel('x$_{Dot}$ / $\pi$')
ax.grid()
#ax.set_ylim([0,1.1])
ax.plot(X/pi,Xdot/pi,'b',lw = 2)
ax.plot(X[k[0]]/pi,0,'bo',ms = 6)
ax.plot(X[k[1]]/pi,0,'ro',ms = 6)
# ax.plot(0.47,0,'bo',ms = 8)
# ax.plot(0.51,0,'ro',ms = 8)
fig2.tight_layout()

#%%  w > a   not fixed points, non uniform motion, period T
num = 599
X = linspace(0,2*pi,num)
A = linspace(0.1,0.9999,num)
T = np.zeros(num)
for c in range(num):
    fn = 1 / (w - A[c]*sin(X))
    T[c] = integrate.simpson(fn,X)

# FIGURE 3  a vs T
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('a')
ax.set_ylabel('T')
ax.set_title('$\omega$ = 1 > a')
ax.grid()
ax.plot(A,T,'b',lw = 2)
fig3.tight_layout()

#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')