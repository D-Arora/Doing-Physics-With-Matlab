# -*- coding: utf-8 -*-
'''
ds2201.py      Oct 2025
DYNAMICAL SYSTEMS
GRADIENT SYSTEMS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS.ds2200.pdf

Displacement  x
Velocity      v = xDot = y

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time

tStart = time.time()
plt.close('all')

#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x, y = state
    dx = -(3*x**2 - 12*y)
    dy = -(3*y**2 - 12*x)
    return [dx, dy]  

#%%  INPUTS >>>
# Initial conditions xI yI (vI)
xI,yI = -3.88,9
# time span
tS = 0.52; nT = 99999
# Plot color
col = [0,0,1]

#%% Solution ODE for x and y 
t = linspace(0,tS,nT)
u0 = [xI,yI]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]       

#%% PHASE PORTRAIT streamplot
N = 11
x = linspace(-10,10,N)
xx,yy = np.meshgrid(x,x)
xxDot = -(3*xx**2 - 12*yy)
yyDot = -(3*yy**2 - 12*xx)

#%% NULLCLINES
xNx = linspace(-9,9,599)
yNx = xNx**2/4
yNy = linspace(-9,9,599)
xNy = yNy**2/4

#%%   FIG 1: t vs x   t vs y = v
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,4)
fig1, ax = plt.subplots(nrows=2, ncols=1)

C = 0   
ax[C].set_xlabel('t'); ax[C].set_ylabel('x')
ax[C].grid()
ax[C].plot(t,xS,lw = 2,color = col)

C = 1   
ax[C].set_xlabel('t'); ax[C].set_ylabel('y')
ax[C].grid()
ax[C].plot(t,yS,lw = 2,color = col)

fig1.tight_layout()
fig1.savefig('a1.png')

#%% FIGURE 2: Phase Portrait  quiver
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,4)
fig2, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x')
axes.set_ylabel('v', rotation = 0)
axes.grid()
axes.set_aspect('equal', 'box')
axes.set_xlim([-10,10]);axes.set_ylim([-10,10]);
axes.quiver(xx,yy,xxDot,yyDot,linewidth = 1)
axes.plot(xS,yS,'k',lw = 2)
axes.plot(0,0,'mo',ms = 8)
axes.plot(xS[0],yS[0],'go',ms = 8)
axes.plot(4,4,'mo',ms = 8)
axes.plot(xNx,yNx,'b')
axes.plot(xNy,yNy,'r')

N = 99
x = linspace(-10,10,N)
xx,yy = np.meshgrid(x,x)
xxDot = -(3*xx**2 - 12*yy)
yyDot = -(3*yy**2 - 12*xx)
V = xx**3 + yy**3 - 12*xx*yy
plt.contour(xx,yy,V,levels=20)

fig2.tight_layout()
fig2.savefig('a2.png')

#%% FIGURE 3: Phase Portrait  streamplot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,5) #(4,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x')
axes.set_ylabel('y', rotation = 0)

axes.grid()
axes.set_aspect('equal', 'box')

plt.contourf(xx,yy,V,levels=20)
axes.streamplot(xx,yy,xxDot,yyDot,linewidth = 1, density = 0.8, color='w')
axes.plot(0,0,'mo',ms = 6)
axes.plot(4,4,'mo',ms = 6)
#plt.contour(xx,yy,V,levels=20)
plt.colorbar()

fig3.tight_layout()
fig3.savefig('a3.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


