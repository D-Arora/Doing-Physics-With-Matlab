# -*- coding: utf-8 -*-

"""
cs100super.py
Aug 25

DYNAMICAL SYSTEMS:
    Saddle Node Bifurcations: supercritical

# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L6.pdf


"""

# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace, zeros

plt.close('all')
tStart = time.time()


#%%  SOLVE ODE
# Solve ODE for x
def lorenz(t, state):    
    xS = state
    dx = r - xS**2
    return dx  

#%% INPUT r value  r = -16, 0, 16 
r = -16
# INPUT initial condition x(0) = x0 
x0 = 4
# INPUT max time span (sensitive to r and x0 values)
tMax = 0.57

#if r > 0: tMax = 0.19

# x, Xdot and fixed points, plot x vs xDot
x = linspace(-8,8,299)
xDot = r - x**2
if r >= 0:
   xss = zeros(2)
   xss[0] = -(r)**0.5; xss[1] = (r)**0.5

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(1)
ax.set_xlabel('x',color = 'black')
ax.set_ylabel('x$_{Dot}$',color= 'black') 
ax.grid()
ax.set_xlim([-8, 8])
ax.set_xticks(np.arange(-8,8.1,2))
ax.set_ylim([-20, 20])
ax.set_yticks(np.arange(-20,21,5))
if r < 0:
   ax.set_ylim([-80, 0])
   ax.set_yticks(np.arange(-80,1,20)) 
   ax.arrow(2,-10,-4,0,lw = 2,color = 'k', width = 0.2, head_width = 3,
       head_length = 0.4)

ax.set_title(r'r = %2.0f' %r,  fontsize = 12)
ax.plot(x,xDot,'b',lw = 2)

if r > 0:
   xP = xss[0]; yP = 0; col = [1,0,0]
   ax.plot(xP,yP,'o',color = col, ms = 8)
   xP = xss[1]; yP = 0; col = [0,0,1]
   ax.plot(xP,yP,'o',color = col, ms = 8)
   ax.arrow(-5,0,-1,0,lw = 2,color = 'k', width = 0.2, head_width = 3,
          head_length = 0.4)
   ax.arrow(-3,0,1,0,lw = 2,color = 'k', width = 0.2, head_width = 3,
          head_length = 0.4)
   ax.arrow(6,0,-1,0,lw = 2,color = 'k', width = 0.2, head_width = 3,
          head_length = 0.4)

if r == 0:
   xP = -0.2; yP = 0; col = [1,0,0]
   ax.plot(xP,yP,'o',color = col, ms = 8)
   xP = 0.2; yP = 0; col = [0,0,1]
   ax.plot(xP,yP,'o',color = col, ms = 8)
   ax.arrow(-4,0,-1,0,lw = 2,color = 'k', width = 0.2, head_width = 3,
          head_length = 0.4)
   ax.arrow(4,0,-1,0,lw = 2,color = 'k', width = 0.2, head_width = 3,
          head_length = 0.4)
  
fig1.tight_layout()


#%% SOLVE ODE
# time span
N = 9999; 
t = linspace(0,tMax,N)

sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3,2.5)
fig2, ax = plt.subplots(1)
ax.set_xlabel('t',color = 'black')
ax.set_ylabel('x',color= 'black') 
ax.set_title(r'r = %2.0f' %r  + '   x(0) = %2.3f' %x0, fontsize = 12)
ax.grid()

# ax.set_ylim([y1, y2])
# ax.set_yticks(np.arange(y1,y2,dy))

xP = t; yP = xS; col = [0,0,1]
ax.plot(xP,yP,'b',lw = 2)    
       
fig2.tight_layout()


#%% BIFURCATION DIAGRAM
R = linspace(0,20,599)            # r values
xe1 = (R)**0.5                    # fixed points
xe2 = -xe1
y1 = -5; y2 = 5; dy = 1

fig3, ax = plt.subplots(1)
ax.set_xlabel('r',color = 'black',fontsize = 12)
ax.set_ylabel('$x_e$',color= 'black',fontsize = 14)
ax.grid()
ax.set_ylim([y1, y2])
ax.set_xlim([0,21])
ax.set_xticks(np.arange(0,20.1,4))
xP = R; yP = xe1; col = [0,0,1]
ax.plot(xP,yP,color = col,lw = 2,label = 'unstable') 
xP = R; yP = xe2; col = [1,0,0]
ax.plot(xP,yP,color = col,lw = 2,label = 'stable') 
xP = 16; yP = -4; col = [1,0,0]
ax.plot(xP,yP,'o',color = col,ms = 8) 
xP = 16; yP = 4; col = [0,0,1]
ax.plot(xP,yP,'o',color = col,ms = 8) 
ax.legend(fontsize = 10) 
fig3.tight_layout()
 

#%% 
fig1.savefig('a1.png')     
fig2.savefig('a2.png')
fig3.savefig('a3.png')