# -*- coding: utf-8 -*-
'''
# ds1501.py     Sep 2025
# NONLINEAR [2D] DYNAMICAL SYSTEMS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS/ds1501.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, real, imag, exp 
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
    dx = x + exp(-y)
    dy = -y
    return [dx, dy]  

#%% INPUTS
# X and Y limits for plots 
L1 = -6; L2 = 6; nL = 599
# Time span 
t1 = 0; nT = 599

# Initial values: a, x0, y0    comment / uncomment for different inputs
N = 9    # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

x0[0], y0[0] = 1.00, 4.50
x0[1], y0[1] = -2.90, 1.31 
x0[2], y0[2] = -3.00, -4.00
x0[3], y0[3] = 1.00, -2.00
x0[4], y0[4] = 3.00, -1.50
x0[5], y0[5] = 1.50, -0.50
x0[6], y0[6] = 0.10, 0.10
t2 = 5


#%% SETUP

# Solution ODE for x and y 
num = len(x0)
xS = zeros([nT,num]) ; yS = zeros([nT,num])
t = linspace(t1,t2,nT)

for c in range(num):
  u0 = [x0[c],y0[c]]
  sol = odeint(lorenz, u0, t, tfirst=True)
  xS[:,c] = sol[:,0]     
  yS[:,c] = sol[:,1] 
  
# Phase Portrait: quiver and stream plots  
xP = linspace(L1,L2,8); yP = xP
xx,yy = np.meshgrid(xP,yP)
xxDot = xx + exp(-yy)
yyDot = -yy

# Nullclines
yR = linspace(L1,L2,nL)
xR = linspace(L1,L2,nL)
xN = -exp(-yR) 
yN = zeros(nL)

dfdx = 1; dfdy = -1
dgdx = 0; dgdy = -1
J = [[dfdx,dfdy],[dgdx,dgdy]]
eV, eF = eig(J)


# Manifolds
# xM0 = JeigF[0,0]; yM0 = JeigF[1,0]
# xM1 = JeigF[0,1]; yM1 = JeigF[1,1]
# m0 = yM0/xM0   
# m1 = yM1/xM1
# XM = xP
# YM0 = m0*XM; YM1 = m1*XM

    
#%% Console output
# print(' ')
# print('A matrix: a11 = %0.2f' %a11 + '  a12 = %0.2f' %a12
#       + '  a11 = %0.2f ' %  a21 + '  a12 = %0.2f' %a22)
# print('Determinant A = %0.5f' % detA)
# print('Initial conditions (x0, y0)')
# for c in range(num):
#     print('  (%0.4f, ' %x0[c] + '%0.4f)' %y0[c])
# print('Eigenvalues    ', JeigV[0],'   ', JeigV[1])
# print('Eigenfunctions ', JeigF[0],'   ', JeigF[1]) 


#%% GRAPHICS
# FIGURE 1: x vs xDot   y vs yDot
plt.rcParams['font.size'] = 11
plt.rcParams["figure.figsize"] = (7,3)
fig1, ax = plt.subplots(nrows=1, ncols=2)
XP = linspace(L1,L2,88); YP = XP
XX,YY = np.meshgrid(XP,XP)
XXDot = XX + exp(-YY)
YYDot = -YY

C = 0   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('x$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([L1,L2])
cf = ax[C].pcolor(XX,YY,XXDot, cmap='jet') 
fig1.colorbar(cf, ax=ax[C])

C = 1   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('y$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([L1,L2])
ax[C].xaxis.grid()
ax[C].yaxis.grid()
cf = ax[C].pcolor(XX,YY,YYDot, cmap='jet') 
fig1.colorbar(cf, ax=ax[C])
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

C = 1   
axes[C].set_xlabel('t',color= 'black',fontsize = 12)
axes[C].set_ylabel('y',color = 'black',fontsize = 12)
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(t, yS,lw = 2)
fig2.tight_layout()



#%% FIGURE 3: Phase Portrait  quiver plot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
axes.xaxis.grid()
axes.yaxis.grid()
axes.set_xlim([L1,L2])
axes.set_ylim([L1,L2])
K = np.sqrt(xxDot**2 + yyDot**2)
axes.quiver(xx,yy,xxDot/K,yyDot/K)
axes.plot(xS,yS,'g',lw = 2)

axes.plot(xN,yR,'b',lw = 2)
axes.plot(xR,yN,'r',lw = 2)

axes.plot(x0,y0,'go', ms = 6)
axes.plot(-1,0,'ko',ms = 6)

fig3.tight_layout()

#%% FIGURE 4: Phase Portrait  streamine
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig4, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
axes.xaxis.grid()
axes.yaxis.grid()
axes.set_xlim([L1,L2])
axes.set_ylim([L1,L2])
axes.streamplot(xx,yy,xxDot,yyDot)
axes.plot(xS,yS,'g',lw = 2)

axes.plot(xN,yR,'b',lw = 2)
axes.plot(xR,yN,'r',lw = 2)

axes.plot(x0,y0,'go', ms = 6)
axes.plot(-1,0,'ko',ms = 6)

fig4.tight_layout()

# #%%  FIGURE 5: Phase Portrait  quiver plot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,3)
fig5, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
axes.xaxis.grid()
axes.yaxis.grid()
axes.set_xlim([L1,L2])
axes.set_ylim([L1,L2])

# K = np.sqrt(xxDot**2 + yyDot**2)
# axes.quiver(xx,yy,xxDot/K,yyDot/K)
# axes.plot(xS,yS,'g',lw = 2)

axes.plot(xN,yR,'b',lw = 2,label = 'x-null')
axes.plot(xR,yN,'r',lw = 2,label = 'y-null')

axes.arrow(3, 0, 1, 0, head_width=0.4, head_length=0.5, lw = 3,fc='k', ec='k')
axes.arrow(-2, 0, -1, 0, head_width=0.4, head_length=0.5, lw = 3,fc='k', ec='k')
axes.plot(-1,0,'ko',ms = 6)

yA = -1.5; xA = -exp(-yA)
axes.arrow(xA, yA, 0, 1, head_width=0.4, head_length=0.5, lw = 3,fc='k', ec='k')

yA = 4; xA = -exp(-yA)
dy = 1; dx = exp(-yA)
axes.arrow(xA, yA, 0, 1, head_width=0.4, head_length=0.5, lw = 3,fc='k', ec='k')


axes.set_aspect('equal', adjustable='box')
axes.legend()
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
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')



'''
#