# -*- coding: utf-8 -*-
'''
ds1800.py      Oct 2025
[2D] NON-LINEAR DYNAMICAL SYSTEMS
CONSERVATIVE SYSTEMS: DOULBLE POTENTIAL WELL

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS1800.pdf

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
    dy = x - x**3
    return [dx, dy]  

def PE(x):
    E = -0.5*x**2 + 0.25*x**4
    return E

#%% POTENTIAL ENERGY FUNCTION
N,L = 999, 1.6; x = linspace(-L,L,N)
V = PE(x)

# FIXED POINTS xE yE
xE = array([-1,0,1]); yE = array([0,0,0])

#%% Solution ODE for x and y 
x0,y0 = 1.1,0.70
t1 = 0; t2 = 30; nT = 9999
t = linspace(t1,t2,nT)
u0 = [x0,y0]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]       

KS = 0.5*yS**2
VS = PE(xS)
ES = KS+VS


#%% PHASE PORTRAIT streamplot
xx,yy = np.meshgrid(x,x)
xxDot = yy
yyDot = xx - xx**3


#%% GRAPHICS
# FIGURE 1: potetnail energy function V(x)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6.5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x'); ax.set_ylabel('V    E')
ax.grid()
ax.set_xlim([-1.6,1.6])
ax.plot(x,V,'b',lw = 2)
yP = 0.02*ones(N); ax.plot(x,yP,'k',lw = 1)
yP = -0.02*ones(N); ax.plot(x,yP,'k',lw = 1)
yP = -0.22*ones(N); ax.plot(x,yP,'k',lw = 1)
yP = 0.2*ones(N); ax.plot(x,yP,'k',lw = 1)
fig1.tight_layout()
fig1.savefig('a1.png')


# FIGURE 2: t vs x   t vs y = v
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,6)
fig2, ax = plt.subplots(nrows=3, ncols=1)
fig2.subplots_adjust(top=0.97, bottom = 0.10, left = 0.15,\
                  right = 0.97, hspace = 0.75,wspace=0.2)

C = 0   
ax[C].set_xlabel('t')
ax[C].set_ylabel('x')
ax[C].grid()
ax[C].plot(t,xS,'b',lw = 2)

C = 1   
ax[C].set_xlabel('t')
ax[C].set_ylabel('v')
ax[C].grid()
ax[C].plot(t,yS,'r',lw = 2)

C = 2   
ax[C].set_xlabel('t')
ax[C].set_ylabel('K  V  E')
ax[C].set_title('E = %0.5f' % max(ES))
ax[C].grid()
ax[C].plot(t,KS,'r',lw = 2)
ax[C].plot(t,VS,'b',lw = 2)
ax[C].plot(t,ES,'k',lw = 2)

fig2.tight_layout()
fig2.savefig('a2.png')


#%% FIGURE 3: Phase Portrait  streamplot
# Phase Portrait: quiver and stream plots  
plt.rcParams['font.size'] = 10
#plt.rcParams["figure.figsize"] = (7,5)
plt.rcParams["figure.figsize"] = (7,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x')
axes.set_ylabel('v')
axes.grid()

axes.streamplot(xx,yy,xxDot,yyDot,linewidth = 1, density = 1.2)
axes.plot(xE,yE,'ko',ms = 6)

fig3.tight_layout()
fig3.savefig('a3.png')

#xxx
#%% SOLVE ODE
# Time span 
t1 = 0; t2 = 40; nT = 599
t = linspace(t1,t2,nT)
N = 7    # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

x0[0], y0[0] = 1.1,1.1
x0[1], y0[1] =1.2,0
x0[2], y0[2] = 1.2,0.7
x0[3], y0[3] = 1.2,0.6
x0[4], y0[4] = -1,-0.65
x0[5], y0[5] = -1,0.1
x0[6], y0[6] = -1.8,0.1

xS = zeros([nT,N]); yS = xS

'''
x0 = 1.1; y0 = 0.6913392799487094
u0 = 1.1, 0.6913392799487094
'''

# for c in range(N):
#   u0 = [x0[c],y0[c]]
#   sol = odeint(lorenz, u0, t, tfirst=True)
#   xS = sol[:,0]     
#   yS = sol[:,1]       
#   plt.plot(xS,yS,lw = 2.5)
#   plt.plot(xS[0],yS[0],'go',ms = 6) 


# y0 = 0.6913392799487094
# V = PE(x0)
# K = 0.5*y0**2
# Etot = K+V

# print('   x0    y0      E')
# for c in range(N):
#     print('  %0.2f' %x0[c] + '   %0.2f' %y0[c] + '   %0.2f  ' %Etot[c] )
    
    
fig3.savefig('a3A.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


