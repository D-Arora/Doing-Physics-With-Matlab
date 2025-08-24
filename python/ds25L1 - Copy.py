# -*- coding: utf-8 -*-
"""
ds25L.py
Aug 25

DYNAMICAL SYSTEMS:
    SIMPLE AND DAMPED HARMONICS MOTION
    MASS-SPRING SYSTEM


# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L1.pdf

"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace, sqrt

plt.close('all')

tStart = time.time()


#%%  SOLVE ODE x
def lorenz(t, state):    
    x1, x2 = state
    dx1 = x2
    dx2 = -(b/m)*x2 - (k/m)*x1 + AD*cos(wD*t)
    return dx1, dx2 
#%% SETUP
u0 = [0.1,0]
tMax = 10; N = 9999
m = 2
b = 0
k = 36
w = sqrt(k/m)
T = 2*pi/w
f = 1/T

fD = 0.69
AD = 1
wD = 2*pi*fD
TD = 1/fD

#%% SOLVE ODE
t = linspace(0,tMax,N)
sol = odeint(lorenz, u0, t, tfirst=True)
x = sol[:,0] 
v = sol[:,1]
a = -(b/m)*v - (k/m)*x
F = m*a


#%% FIGURE 1: t vs x
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (8,8)
fig1, ax = plt.subplots(nrows=3, ncols=2)
plt.suptitle('m = %0.1f kg' %m + '   b = %0.1f kg/s' %b + '   k = %0.1f N/m' %k +
             '\n   T$_0$ = %0.2f s' %T + '  f$_0$ = %0.2f Hz ' %f
             + '   T$_D$ = %0.2f s' %TD + '   A$_D$ = %0.2f m/s' % AD)
R = 0; C = 0
ax[R,C].set_xlabel('t')
ax[R,C].set_ylabel('x')
ax[R,C].grid()
ax[R,C].plot(t,x,'b',lw = 2)

R = 0; C = 1
ax[R,C].set_xlabel('t')
ax[R,C].set_ylabel('v')
ax[R,C].grid()
ax[R,C].plot(t,v,'b',lw = 2)

R = 1; C = 0
ax[R,C].set_xlabel('t')
ax[R,C].set_ylabel('a')
ax[R,C].grid()
ax[R,C].plot(t,a,'b',lw = 2)

R = 1; C = 1
ax[R,C].set_xlabel('t')
ax[R,C].set_ylabel('F')
ax[R,C].grid()
ax[R,C].plot(t,F,'b',lw = 2)

R = 2; C = 0
ax[R,C].set_xlabel('x')
ax[R,C].set_ylabel('F')
ax[R,C].grid()
ax[R,C].plot(x,F,'b',lw = 2)

R = 2; C = 1
ax[R,C].set_xlabel('x')
ax[R,C].set_ylabel('v')
ax[R,C].grid()
ax[R,C].plot(x,v,'b',lw = 2)

fig1.tight_layout()




#%%
fig1.savefig('a1.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
