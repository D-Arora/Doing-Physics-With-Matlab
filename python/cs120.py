# -*- coding: utf-8 -*-
'''
# cs120.py      Sep 2025
# NONLINEAR [2D] DYNAMICAL SYSTEMS
# SADDLE NODE BIFURCATION r = 9

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS1503.pdf

'''

# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, sqrt, ones
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')

tStart = time.time()

#%% FUNCTIONS  Solve ODE for x,y    x = theta   y = omega 
def lorenz(t, state):    
    x, y = state
    dx = r - x**2
    dy = -y
    return [dx, dy]  

r = 9                 # control parameter
u0 = [-2.99,5]       # initial conditions 
t2 = 1               # simulation time
L = 6; nL = 599       # system dimensions

#%% x vs xDot   y vs yDot 
x = linspace(-L,L,nL); y = linspace(-L,L,nL)
xDot = r - x**2; yDot = - y

# Solution ODE for x and y  
t1 = 0; nT = 999
t = linspace(t1,t2,nT)
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]            

# Phase Portrait quiver plot  
x1 = -6; x2 = 6; nX = 9
xP = linspace(x1,x2,nX)
y1 = -6; y2 = 6; nY = 9
yP = linspace(y1,y2,nX)
xx,yy = np.meshgrid(xP,yP)
xxDot = r - xx**2
yyDot = -yy

# Jacobian matrix and eigenvalues
J = np.array([[0,0],[0,-1]])
Jev, Jef = eig(J)
print(Jev)

# fixed points
xe = zeros(2); ye = zeros(2)
xe[0], ye[0] = -sqrt(r), 0
xe[1], ye[1] = sqrt(r), 0

#%% GRAPHICS
# FIGURE 1: x vs xDot   y vs yDot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,3)
fig1, axes = plt.subplots(nrows=1, ncols=2)
    
C = 0   
axes[C].set_xlabel('x',color= 'black',fontsize = 12)
axes[C].set_ylabel('$x_{dot}$',color = 'black',fontsize = 12)
axes[C].set_title('r = %2.1f' % r, fontsize = 14)
axes[C].set_xlim([-5, 5])
axes[C].set_xticks(np.arange(-5,5.1,1))
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(x, xDot,'b',lw = 2)
axes[C].plot(-np.sqrt(r),0,'ro',ms = 6)
axes[C].plot(np.sqrt(r),0,'bo',ms = 6)

C = 1   
axes[C].set_xlabel('y',color= 'black',fontsize = 12)
axes[C].set_ylabel('$y_{dot}$',color = 'black',fontsize = 12)
axes[C].set_title('r = %2.1f' % r, fontsize = 14)
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(y, yDot,'b',lw = 2)
axes[C].plot(0,0,'bo',ms = 6)

fig1.tight_layout()


#%% FIGURE 2: t vs x   t vs y
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,3)
fig2, axes = plt.subplots(nrows=1, ncols=2)
    
C = 0   
axes[C].set_xlabel('t',color= 'black',fontsize = 12)
axes[C].set_ylabel('x',color = 'black',fontsize = 12)
axes[C].set_title('r = %2.1f' % r, fontsize = 14)
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(t, xS,'b',lw = 2)

C = 1   
axes[C].set_xlabel('t',color= 'black',fontsize = 12)
axes[C].set_ylabel('y',color = 'black',fontsize = 12)
axes[C].set_title('r = %2.1f' % r, fontsize = 14)
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(t, yS,'b',lw = 2)
fig2.tight_layout()


#%% FIGURE 3: xe vs r   fixed points as a function of r
xP = linspace(0,16,599)
yPp = np.sqrt(xP); yPm = -yPp

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('r',color= 'black',fontsize = 12)
axes.set_ylabel('$x_e$',color = 'black',fontsize = 12)
axes.xaxis.grid()
axes.yaxis.grid()
axes.plot(xP, yPp,'b',lw = 2)
axes.plot(xP, yPm,'r',lw = 2)
fig3.tight_layout()


#%% FIGURE 4: Phase Portrait  streamplot plot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig4, axes = plt.subplots(nrows=1, ncols=1)

axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)
axes.set_xlim([-6.1, 6.1])
axes.xaxis.grid()
axes.yaxis.grid()

# Nullclines
xP = x; yP = zeros(nL)
axes.plot(xP,yP,'m')
xP = ones(nL)*xe[0]; yP = y
axes.plot(xP,yP,'g')
axes.plot(-xP,yP,'g')
# Fixed points
axes.plot(xe[0],ye[0],'ro',ms = 8)
axes.plot(xe[1],ye[1],'bo',ms = 8)
# Vector field
axes.streamplot(xx,yy,xxDot,yyDot)
fig4.tight_layout()


#%% FIGURE 5: Phase Portrait  quiver plot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig5, axes = plt.subplots(nrows=1, ncols=1)

  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)
axes.set_xlim([-L, L])
axes.xaxis.grid()
axes.yaxis.grid()
# Nullclines
xP = x; yP = zeros(nL)
axes.plot(xP,yP,'m')
xP = ones(nL)*xe[0]; yP = y
axes.plot(xP,yP,'g')
axes.plot(-xP,yP,'g')
# Fixed points
axes.plot(xe[0],ye[0],'ro',ms = 8)
axes.plot(xe[1],ye[1],'bo',ms = 8)
# Vector field
axes.quiver(xx,yy,xxDot,yyDot)
fig5.tight_layout()


#%%  FIGURE 6  Phase portrait + trajectories
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig6, axes = plt.subplots(nrows=1, ncols=1)

axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)
axes.set_xlim([-6.1, 6.1])
axes.xaxis.grid()
axes.yaxis.grid()
# Fixed points
axes.plot(xe[0],ye[0],'ro',ms = 8)
axes.plot(xe[1],ye[1],'bo',ms = 8)
# Vector field
axes.streamplot(xx,yy,xxDot,yyDot)
# Trajectories
x0,y0 = 5,-5
t2 = 2
t = linspace(t1,t2,nT)
u0 = [x0,y0]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]    
axes.plot(xS,yS,'k',lw = 3)
fig6.tight_layout()


#%%   FIGURE 7: slope angle dy/dx
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)
fig7, ax = plt.subplots(nrows=1, ncols=1)
N= 555
XP = linspace(-L,L,N); YP = XP
XX,YY = np.meshgrid(XP,XP)
dX = r - XX**2
dY = -YY
theta = np.arctan2(dY, dX)
ZZ = theta / pi
ax.set_xlabel('x',color= 'black',fontsize = 12)
ax.set_ylabel('y',color = 'black',fontsize = 12)
ax.set_title('slope angle [ rad / $\pi$ ]',color = 'black',fontsize = 14)
ax.set_xticks([-6,-3,0,3,6])
ax.set_yticks([-6,-3,0,3,6])
cf = ax.pcolor(XX,YY,ZZ, cmap='jet') 
# Nullclines
xP = x; yP = zeros(nL)
ax.plot(xP,yP,'w', lw = 1)
xP = ones(nL)*xe[0]; yP = y
ax.plot(xP,yP,'w',lw=1)
ax.plot(-xP,yP,'w',lw=1)
ax.set_box_aspect(1)
cbar = plt.colorbar(cf, ax=ax)
cbar.set_ticks([-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0])
fig7.tight_layout()


#%% IGURE 8: x vs xDot   y vs yDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,4)
fig8, ax = plt.subplots(nrows=1, ncols=2)

C = 0   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('x$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([-L,L])
cf = ax[C].pcolor(XX,YY,dX, cmap='jet') 
# Nullclines
xP = x; yP = zeros(nL)
ax[C].plot(xP,yP,'w', lw = 1)
xP = ones(nL)*xe[0]; yP = y
ax[C].plot(xP,yP,'w',lw=1)
ax[C].plot(-xP,yP,'w',lw=1)
ax[C].set_xticks([-6,-3,0,3,6])
ax[C].set_yticks([-6,-3,0,3,6])
ax[C].set_box_aspect(1)
fig8.colorbar(cf, ax=ax[C])

C = 1   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('y$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([-L,L])
cf = ax[C].pcolor(XX,YY,dY, cmap='jet') 
fig8.colorbar(cf, ax=ax[C])
# Nullclines
xP = x; yP = zeros(nL)
ax[C].plot(xP,yP,'w', lw = 1)
xP = ones(nL)*xe[0]; yP = y
ax[C].plot(xP,yP,'w',lw=1)
ax[C].plot(-xP,yP,'w',lw=1)
ax[C].set_xticks([-6,-3,0,3,6])
ax[C].set_yticks([-6,-3,0,3,6])
ax[C].set_box_aspect(1)
fig8.tight_layout()


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

#%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')
fig6.savefig('a6.png')
fig7.savefig('a7.png')
fig8.savefig('a8.png')

'''
