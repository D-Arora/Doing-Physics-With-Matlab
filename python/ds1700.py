# -*- coding: utf-8 -*-
'''
ds1700.py      Oct 2025
[2D] NON-LINEAR DYNAMICAL SYSTEMS
RABBITS and SHEEP population dynamics

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS1700.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, real, imag 
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
    dx = x*(3-x-2*y)
    dy = y*(2-x-y)
    return [dx, dy]  


#%% FIXED POINTS xE yE
xE = array([0,0,3,1]); yE = array([0,2,0,1])

def JM(x,y):
    a = zeros([2,2])
    a[0,0] = 3 - 2*x - 2*y
    a[0,1] = -2*x
    a[1,0] = -y
    a[1,1] = 2 - x - 2*y
    return a
# Jacobian
J00 = JM(xE[0],yE[0])
J02 = JM(xE[1],yE[1])
J30 = JM(xE[2],yE[2])
J11 = JM(xE[3],yE[3])
#eigenvalues v  and eigenvectors V

v00, V00 = eig(J00)
v02, V02 = eig(J02)
v30, V30 = eig(J30)
v11, V11 = eig(J11)

# NULLCLINES
num = 599
X = linspace(0,3,num)
xN = (3 - X)/2
yN = 2 - X

print('Eigenvalues v')
print('v00 = ', v00)
print('v02 = ', v02)
print('v30 = ', v30)
print('v11 = ', v11)
print('\n  Eigenvectors V')
print('V00 = ', V00)
print('V02 = ', V02)
print('V30 = ', np.round(V30,2))
print('V11 = ', V11)


#%% Solution ODE for x and y 
x0,y0 = 3,1
t1 = 0; t2 = 30; nT = 9999
t = linspace(t1,t2,nT)
u0 = [x0,y0]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]       


#%% GRAPHICS
# FIGURE 1: x vs xDot   y vs yDot
L = 4  # Phase space dimensions
plt.rcParams['font.size'] = 11
plt.rcParams["figure.figsize"] = (7,3)
fig1, ax = plt.subplots(nrows=1, ncols=2)
XP = linspace(-0.10,L,88); YP = XP
XX,YY = np.meshgrid(XP,XP)
XXDot = XX*(3-XX-2*YY)
YYDot = YY*(2-XX-YY)
C = 0   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('x$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([0,L]); ax[C].set_ylim([0,L])
cf = ax[C].pcolor(XX,YY,XXDot, cmap='hot') 
ax[C].plot(X,xN,'b'); ax[C].plot(X,yN,'m')
ax[C].plot(xE,yE,'ko',ms = 8)
fig1.colorbar(cf, ax=ax[C])

C = 1   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('y$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([0,L]); ax[C].set_ylim([0,L])
cf = ax[C].pcolor(XX,YY,YYDot, cmap='hot') 
ax[C].plot(X,xN,'b'); ax[C].plot(X,yN,'m')
ax[C].plot(xE,yE,'ko',ms = 8)
fig1.colorbar(cf, ax=ax[C])
fig1.tight_layout()
fig1.savefig('a1.png')


#%% FIGURE 2: Phase Portrait  streamplot
# Phase Portrait: quiver and stream plots  
xP = linspace(0,4,99); yP = xP
xx,yy = np.meshgrid(xP,yP)
xxDot = xx*(3-xx-2*yy)
yyDot = yy*(2-xx-yy)
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,4)
fig2, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
#axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
axes.grid()
axes.set_xlim([0,4])
axes.set_ylim([0,4])
axes.streamplot(xx,yy,xxDot,yyDot)
axes.plot(xE,yE,'ko',ms = 6)

axes.set_aspect('equal', 'box')
fig2.tight_layout()


#%% SOLVE ODE
# Time span 
t1 = 0; t2 = 510; nT = 599
t = linspace(t1,t2,nT)
N = 6    # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

x0[0], y0[0] = 3.6,2.5
x0[1], y0[1] = 3.6,2.6
x0[2], y0[2] = 0.1,0.3
x0[3], y0[3] = 0.1,0.1
x0[4], y0[4] = 3,3.8
x0[5], y0[5] = 3.8,1.5
xS = zeros([nT,N]); yS = xS

for c in range(N):
  u0 = [x0[c],y0[c]]
  sol = odeint(lorenz, u0, t, tfirst=True)
  xS = sol[:,0]     
  yS = sol[:,1]       
  plt.plot(xS,yS,lw = 2.5)
  plt.plot(xS[0],yS[0],'go',ms = 5)

fig2.savefig('a2.png')
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


