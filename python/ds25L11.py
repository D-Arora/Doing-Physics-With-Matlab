# -*- coding: utf-8 -*-
"""
ds25L11.py
Sep 25

DYNAMICAL SYSTEMS:
    Flow on a circle


# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L11.pdf
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
w = 0.8
a = 1
wa = w/a
x = linspace(0,2*pi,num)
xDot = XDOT(x)
A = linspace(0.1,0.9999,num)
T = np.zeros(num)
for c in range(num):
    fn = 1 / (w - A[c]*sin(x))
    T[c] = integrate.simpson(fn,x)

# a vs T
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('a')
ax.set_ylabel('T')
ax.set_title('$\omega$ = 1 > a ')
ax.grid()
ax.plot(A,T,'b',lw = 2)
fig1.tight_layout()
fig1.savefig('a1.png')

# x vs xDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.5)
fig1A, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x / $\pi$')
ax.set_ylabel('x$_{Dot}$ / $\pi$')
ax.set_title('$\omega$ = %0.3f' %w + '   a = %0.3f' %a 
             + '   $\omega$/a = %0.3f' %wa , fontsize = 12)
ax.grid()
ax.set_xlim([0,2])
ax.set_xticks(np.arange(0,2.2,0.25))
ax.plot(x/pi,xDot/pi,'k',lw = 2)
ax.plot([0,2],[0,0],'m',lw = 2)
fig1A.tight_layout()
fig1A.savefig('a1A.png')

# t vs x
N = 9999
tMax = 20
x0 = 0.5
t = linspace(0,tMax,N)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0]
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.5)
fig1T, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_title('$\omega$ = %0.3f' %w + '   a = %0.3f' %a 
             + '   $\omega$/a = %0.3f' %wa , fontsize = 12)
ax.grid()

ax.plot(t,xS/pi,'k',lw = 2)

#ax.set_xlim([0,2])
#ax.set_xticks(np.arange(0,2.2,0.25))

fig1T.tight_layout()
fig1T.savefig('a1T.png')


#%% CELL 2:  w = a
w = 1
a = 1
N = 999
wa = w/a
x = linspace(0,2*pi,N)
xDot = XDOT(x)

# zero xDot
xZ = pi/2

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.5)
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

ax.plot(0.51,0,'ro',ms = 6)
ax.plot(0.49,0,'bo',ms = 6)
fig2.tight_layout()
fig2.savefig('a2.png')

# SOLVE ODE for trajectories
N = 9999
tMax = 20
u0 = np.array([0,1,0.25,0.60])*pi
t = linspace(0,tMax,N)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.5)
fig2T, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_title('$\omega$ = %0.3f' %w + '   a = %0.3f' %a 
             + '   $\omega$/a = %0.3f' %wa , fontsize = 12)
ax.grid()

for c in range(4):
    x0 = u0[c]
    sol = odeint(lorenz, x0, t, tfirst=True)
    xS = sol[:,0]
    ax.plot(t,xS/pi,'k',lw = 2)

#ax.set_xlim([0,2])
#ax.set_xticks(np.arange(0,2.2,0.25))

ax.plot([0,tMax],[0.5,0.5],'b',lw = 2)
ax.plot([0,tMax],[2.5,2.5],'b',lw = 2)
fig2T.tight_layout()
fig2T.savefig('a2T.png')



#%%  CELL 3: w < a
w = 1
a = sqrt(2)
N = 999
wa = w/a
x = linspace(0,2*pi,N)
xDot = XDOT(x)

# Find zeros xDot
Q = np.zeros(2); p = 0
for c in range(N-2):
    q = xDot[c]*xDot[c+1]
    if q <= 0:
       Q[p] = c
       p = int(p+1)     
QI = Q.astype(int)
xZ = x[QI]          # Zeros for xDot  
print(np.round(xZ/pi,3))
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x / $\pi$')
ax.set_ylabel('x$_{Dot}$ / $\pi$')
ax.set_title('$\omega$ = %0.3f' %w + '   a = %0.3f' %a 
             + '   $\omega$/a = %0.3f' %wa , fontsize = 12)
ax.grid()
ax.set_xlim([0,2])
ax.set_xticks(np.arange(0,2.2,0.25))
ax.plot(x/pi,xDot/pi,'b',lw = 2)
ax.plot([0,2],[0,0],'m',lw = 2)
ax.plot(xZ[0]/pi,0,'bo',ms = 6)
ax.plot(xZ[1]/pi,0,'ro',ms = 6)
fig3.tight_layout()
fig3.savefig('a3.png')

# SOLVE ODE for trajectories
N = 9999
tMax = 10
u0 = np.array([0.0,0.4,0.6,0.8,1.0,1.2])*pi
t = linspace(0,tMax,N)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.5)
fig3T, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_title('$\omega$ = %0.3f' %w + '   a = %0.3f' %a 
             + '   $\omega$/a = %0.3f' %wa , fontsize = 12)
ax.grid()

for c in range(6):
    x0 = u0[c]
    sol = odeint(lorenz, x0, t, tfirst=True)
    xS = sol[:,0]
    ax.plot(t,xS/pi,'k',lw = 2)

#ax.set_xlim([0,2])
#ax.set_xticks(np.arange(0,2.2,0.25))

ax.plot([0,tMax],[0.25,0.25],'b',lw = 2)
ax.plot([0,tMax],[0.75,0.75],'r',lw = 2)
ax.plot([0,tMax],[2.25,2.25],'b',lw = 2)
fig3T.tight_layout()
fig3T.savefig('a3T.png')


