# -*- coding: utf-8 -*-
'''
ds1900.py      Oct 2025
[2D] NON-LINEAR DYNAMICAL SYSTEMS
CONSERVATIVE SYSTEMS: SIMPLE DAMPED PENDULUM

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS1900.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import random

tStart = time.time()
plt.close('all')

#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x, y = state
    dx = y
    dy = -b*y - w**2*sin(x)
    return [dx, dy]  

def PE(x):
    E = m*g*L*(1-cos(xS))
    return E

#%%  INPUTS and : Model parameters
m = 1
g = 10
b = 0
T = 1
w = 2*pi/T
L = g/w**2
N = 999
x = linspace(-pi,pi,N)
V = m*g*L*(1 - cos(x))

#%%
omega0 = np.sqrt( (4*g - 2*g*(1-cos(pi/8)))/L )
print(omega0)

# FIXED POINTS xE yE
xE = array([-pi,0,pi]); yE = array([0,0,0])

#%% Solution ODE for x and y 
x0,y0 = pi/8,10

t1 = 0; t2 = 5; nT = 9999
t = linspace(t1,t2,nT)
u0 = [x0,y0]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]       

KS = 0.5*(L*yS)**2
VS = PE(xS)
ES = KS+VS

print(ES[0])

#%% PHASE PORTRAIT streamplot
x = linspace(-1.1*pi,1.1*pi,N)
xx,yy = np.meshgrid(x,4*x)
xxDot = yy
yyDot = -b*yy - w**2*sin(xx)

#%%  FIG 1: potential energy function V(x)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6.5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel(r'$\theta  /  \pi$'); ax.set_ylabel('V  [ J ]')
ax.grid()
ax.set_xlim([-1,1])
ax.plot(x/pi,V,'b',lw = 2)

fig1.tight_layout()
fig1.savefig('a1.png')

#%%   FIG 2: t vs x   t vs y = v
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,6)
fig2, ax = plt.subplots(nrows=3, ncols=1)
fig2.subplots_adjust(top=0.97, bottom = 0.10, left = 0.15,\
                  right = 0.97, hspace = 0.75,wspace=0.2)

C = 0   
ax[C].set_xlabel('t  [ s ]')
ax[C].set_ylabel(r'$\theta / \pi$')
ax[C].grid()
ax[C].plot(t,xS/pi,'b',lw = 2)

C = 1   
ax[C].set_xlabel('t  [ s ]')
ax[C].set_ylabel(r'$\omega$  [ rad/s ]')
ax[C].grid()
ax[C].plot(t,yS,'r',lw = 2)

C = 2   
ax[C].set_xlabel('t [ s ]')
ax[C].set_ylabel('K  V  E  [ J ]')
# ax[C].set_title('E = %0.5f' % max(ES))
ax[C].grid()
ax[C].plot(t,KS,'r',lw = 2)
ax[C].plot(t,VS,'b',lw = 2)
ax[C].plot(t,ES,'k',lw = 2)

#fig2.tight_layout()
fig2.savefig('a2.png')

#%% FIGURE 3: Phase Portrait  streamplot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,4)
#plt.rcParams["figure.figsize"] = (7,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel(r'$\theta / \pi$')
axes.set_ylabel(r'$\omega$', rotation = 0)

axes.grid()

axes.streamplot(xx/pi,yy,xxDot,yyDot,linewidth = 1, density = 1.2)
axes.plot(xE/pi,yE,'ko',ms = 6)

fig3.tight_layout()
fig3.savefig('a3.png')

#xxx
#%% SOLVE ODE
# Time span 
t1 = 0; t2 = 5; nT = 999
t = linspace(t1,t2,nT)
N = 5    # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

x0[0], y0[0] = -pi/4,3
x0[1], y0[1] = pi/2,2
x0[2], y0[2] = -pi/1.1,0
x0[3], y0[3] = -pi/3,0
x0[4], y0[4] = pi/1.3,4

for c in range(N):
  u0 = [x0[c],y0[c]]
  sol = odeint(lorenz, u0, t, tfirst=True)
  xS = sol[:,0]     
  yS = sol[:,1]       
  plt.plot(xS/pi,yS,lw = 2.5)
  plt.plot(xS[0]/pi,yS[0],'go',ms = 6) 
    
fig3.savefig('a3A.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


