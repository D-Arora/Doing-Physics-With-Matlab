# -*- coding: utf-8 -*-
"""
cs105.py
Aug 25

DYNAMICAL SYSTEMS:
    Pitchfork Bifurcations

# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L8.pdf


"""

# LIBRARIES  ================================================================
import numpy as np
from numpy import pi, sin, cos, linspace
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')
tStart = time.time()

#%%   FUNCTIONS
def funct(r,x):
    xDot = r*x - x**3
    return xDot

#%%  SOLVE ODE for x
def lorenz(t, state):    
    xS = state
    dx = r*xS - xS**3
    return dx  

#%% INPUT r value   r = -16, 0, +16
r = 16

#%% SETUP
u0 = 1.99*pi
tMax = 10; N = 9999

#%% SOLVE ODE
t = linspace(0,tMax,N)
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0] 

N = 999
if r == 0: x1 = -1; x2 = 1
if r < 0: x1 = -4; x2 = 4 
if r > 0: x1 = -6; x2 = 6   
x = linspace(x1,x2,N)
xDot = r*x - x**3

#%%  Fig 1:  Bifurcation diagrams
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('x$_{dot}$')
ax.set_title('r = %0.0f' %r)
ax.grid()
ax.plot(x,xDot,'k',lw = 2)

if r == 0:
   ax.plot(0,0,'bo',ms = 7)
   ax.arrow(-0.85, 0, 0.2, 0,  lw = 2,head_width=0.1, head_length=0.1, fc='k', ec='k')
   ax.arrow(0.85, 0, -0.2, 0, lw = 2,head_width=0.1, head_length=0.1, fc='k', ec='k')
if r == -16:
   ax.plot(0,0,'bo',ms = 7)
   ax.arrow(-3, 0, 1, 0,  lw = 2,head_width=20, head_length=0.5, fc='k', ec='k')
   ax.arrow(3, 0, -1, 0, lw = 2,head_width=20, head_length=0.5, fc='k', ec='k')
if r == 16:
   ax.plot(0,0,'ro',ms = 7) 
   ax.plot(4,0,'bo',ms = 7) 
   ax.plot(-4,0,'bo',ms = 7) 
   ax.arrow(-6, 0, 1, 0,  lw = 2,head_width=20, head_length=0.5, fc='k', ec='k')
   ax.arrow(-2, 0, -1, 0,  lw = 2,head_width=20, head_length=0.5, fc='k', ec='k')
   ax.arrow(2, 0, 1, 0,  lw = 2,head_width=20, head_length=0.5, fc='k', ec='k')
   ax.arrow(6, 0, -1, 0,  lw = 2,head_width=20, head_length=0.5, fc='k', ec='k')
   
fig1.tight_layout()

#%%  Fig 2:  Bifurcation diagrams
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r')
ax.set_ylabel('x$_e$')

ax.grid()

r = linspace(-16,0,N)
xe = np.zeros(N)
ax.plot(r, xe, 'b',lw = 2)

r = linspace(0,16,N)
xe = np.zeros(N)
ax.plot(r, xe, 'r',lw = 2)

r = linspace(0,16,N)
xe = r**0.5
ax.plot(r, xe, 'b',lw = 2)
ax.plot(r, -xe, 'b',lw = 2)

fig2.tight_layout()



#%% TIME EVOLUTION
N = 9999; t1 = 0; t2 = 1
t = linspace(t1,t2,N)
r = 16

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel('x$_{Dot}$',color= 'black',fontsize = 12)
ax.set_xlabel('t',color = 'black',fontsize = 12)
ax.set_title('r = %2.1f' % r, fontsize = 14)
ax.grid()

x0 = -6
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0]  
ax.plot(t, xS, 'b')

x0 = -0.0001
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0]  
ax.plot(t, xS, 'k')

x0 = 0.0001
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0]  
ax.plot(t, xS, 'r')

x0 = 6
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0]  
ax.plot(t, xS, 'm')

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



