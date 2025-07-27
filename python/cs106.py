# -*- coding: utf-8 -*-
"""
NONLINEAR [1D] DYNAMICAL SYSTEMS
FIXED POINTS, STABILITY ANALYSIS, BIFURCATIONS


# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/cs_101.pdf

"""

# LIBRARIES  ================================================================
import numpy as np
from numpy import pi, sin, cos, linspace
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')


#%%  SOLVE ODE for x   where x is the popiulation N
def lorenz(t, state):    
    x = state
    dx = r*x*(1 - x/k)
    return dx  

#%%
r,k = 1.0, 1.0
num, tMax = 999,10
t = linspace(0,tMax,num)

#%%
# FIGURE 1: Membrane potential x vs time t
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('N')
ax.grid()

x0 = np.array([0.01,0.5,2])

for c in range(3):
    u0 = x0[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    xS = sol[:,0] 
    ax.plot(t,xS,'b',lw = 2)

fig1.tight_layout()

#%% FIGURE 2: x vs xDot
x = linspace(0,k,num)
xDot = r*x*(1 - x/k)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x',fontsize = 14)
ax.set_ylabel('x$_{dot}$',fontsize = 14)
ax.grid()
ax.plot(x,xDot,'b',lw = 2)
ax.plot(0,0,'ro',ms=8)
ax.plot(k,0,'bo',ms=8)
fig2.tight_layout()


#%%  FIGURE 3: slope plot  x vs t
N = 15
tQ = linspace(0,10,N)
x = linspace(0,2,N)

f = x*(1-x)

T,X = np.meshgrid(tQ,x)

dX = np.ones([N,N])
F = X*(1-X)
dY = F/(np.sqrt(dX**2 + F**2))
dX = dX/(np.sqrt(dX**2 + F**2))

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_ylim([0,2])
ax.set_xlim([0,10])

ax.quiver(T,X,dX,dY,color = 'b')

for c in range(3):
    u0 = x0[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    xS = sol[:,0] 
    ax.plot(t,xS,'r',lw = 2)

fig3.tight_layout()  

#%%  FIGURE 4: SLOPE FIELD
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_ylim([0,2])
ax.set_xlim([0,10])

ax.streamplot(T,X,dX,dY,color = 'b')

fig4.tight_layout()  

'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')

'''

