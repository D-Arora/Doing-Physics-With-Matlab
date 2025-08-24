# -*- coding: utf-8 -*-
'''
# cs200.py      June 2025
# NONLINEAR [2D] DYNAMICAL SYSTEMS
# FIXED POINTS, STABILITY ANALYSIS, BIFURCATIONS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs200.pdf


https://www.mdpi.com/2079-9292/13/11/2138

'''

#%% Libraries
import numpy as np
from numpy import pi, sqrt, linspace, zeros, array, real, imag
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
    dx = -s*x*y
    dy = s*x*y - r*y
    dx = dx + 0.9*abs(dy)
    return [dx, dy]  

#%% INPUTS
# X and Y limits for plots 
L1 = 0; L2 = 1000; nL = 599

# Time span 
t1 = 0; t2 = 20; nT = 999

# speed at which disease spreads
s = 3e-3
# recovery rate
r = 0.6

# Initial vales x0, y0    comment / uncomment for different inputs
N = 1   # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

# manual IC
x0[0], y0[0] = 1000,2
#x0[1], y0[1] = 750,1
#x0[2], y0[2] = 500,1

P = x0+y0

#%% SETUP
# Solution ODE for x and y 
num = len(x0)
xS = zeros([nT,num]) ; yS = zeros([nT,num]); rS = zeros([nT,num])
t = linspace(t1,t2,nT)

for c in range(num):
  u0 = [x0[c],y0[c]]
  sol = odeint(lorenz, u0, t, tfirst=True)
  xS[:,c] = sol[:,0]     
  yS[:,c] = sol[:,1]
  rS[:,c] = P[c] - xS[:,c] - yS[:,c]
  rS[rS < 0] = 0 
  

  
# Phase Portrait: quiver and stream plots  
xP = linspace(L1,L2,15); yP = xP
xx,yy = np.meshgrid(xP,yP)
xxDot = -s*xx*yy
yyDot = s*xx*yy - r*yy


# Nullclines X [0] / Y [1]
yN0 = zeros(nL); xN0 = linspace(L1,L2,nL) 
yN1 = linspace(L1,L2,nL)
xN1 = (r/s)*np.ones(nL)
    


#%% Console output
# print(' ')
# print('Critical points')
# print('   xC = ', xC,'   yC = ', yC   ) 

# print('Jacobian matrix J0')
# print(J0)
# print('Eigenvalues J0')
# print(np.round(V0,3))     
# print('Eigenvectors J0')
# print(np.round(F0,3))

# print('Jacobian matrix J1')
# print(J1)
# print('Eigenvalues J1')
# print(np.round(V1,3))     
# print('Eigenvectors J1')
# print(np.round(F1,3))

# print('Jacobian matrix J2')
# print(J2)
# print('Eigenvalues J2')
# print(np.round(V1,3))     
# print('Eigenvectors J2')
# print(np.round(F2,3))


#%% GRAPHICS
# FIGURE 1: x vs xDot   y vs yDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,4)
fig1, ax = plt.subplots(nrows=1, ncols=2)
XP = linspace(L1,L2,88); YP = XP
XX,YY = np.meshgrid(XP,XP)
XXDot = -s*XX*YY
YYDot = s*XX*YY - r*YY

C = 0   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('x$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([L1,L2])
cf = ax[C].pcolor(XX,YY,XXDot, cmap='jet') 
#ax[C].plot([L1,L2],[0,0],'w', lw = 2)
ax[C].set_box_aspect(1)
fig1.colorbar(cf, ax=ax[C])

C = 1   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('y$_{dot}$',color = 'black',fontsize = 14)

cf = ax[C].pcolor(XX,YY,YYDot, cmap='jet') 
fig1.colorbar(cf, ax=ax[C])


#ax[C].plot(xN1,yN1,'w', lw = 2)

ax[C].set_box_aspect(1)

ax[C].set_xlim([L1,L2])
ax[C].set_ylim([L1,L2])
ax[C].set_xticks([-2,-1,0,1,2])
ax[C].set_yticks([-2,-1,0,1,2])

fig1.tight_layout()


#%% FIGURE 2: t vs x   t vs y
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6,3)
fig2, axes = plt.subplots(nrows=1, ncols=2)
    
C = 0   
axes[C].set_xlabel('t',color= 'black',fontsize = 12)
axes[C].set_ylabel('x',color = 'black',fontsize = 12)
axes[C].xaxis.grid()
axes[C].yaxis.grid()

axes[C].plot(t, xS,lw = 2)
axes[C].plot(t, rS,lw = 2)

C = 1   
axes[C].set_xlabel('t',color= 'black',fontsize = 12)
axes[C].set_ylabel('y',color = 'black',fontsize = 12)
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(t, yS,lw = 2)
fig2.suptitle('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
fig2.tight_layout()


#%% FIGURE 4: Phase Portrait  quiver plot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,4)
fig4, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
axes.xaxis.grid()
axes.yaxis.grid()
axes.set_xlim([L1,L2])
axes.set_ylim([L1,L2])

axes.quiver(xx,yy,xxDot,yyDot)

axes.plot(xS,yS,'g',lw = 2)
axes.plot(x0,y0,'go', ms = 6)

# axes.plot(xC,yC,'ro', ms = 7)

# axes.plot(xN0,yN0,'r', lw = 1)
# axes.plot(xN1,yN1,'b', lw = 1)

axes.set_box_aspect(1)

fig4.tight_layout()

#%% FIGURE 5: Phase Portrait  streamine
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,4)
fig5, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
axes.xaxis.grid()
axes.yaxis.grid()
axes.set_xlim([L1,L2])
axes.set_ylim([L1,L2])
axes.streamplot(xx,yy,xxDot,yyDot)



# axes.plot(xC,yC,'ro', ms = 7)

# axes.plot(xN0,yN0,'r', lw = 1)
# axes.plot(xN1,yN1,'b', lw = 1)

axes.plot(xS,yS,'g',lw = 3)
axes.plot(x0,y0,'go', ms = 6)

axes.set_box_aspect(1)

fig5.tight_layout()


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


#%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')


'''
#