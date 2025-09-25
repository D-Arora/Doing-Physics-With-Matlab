# -*- coding: utf-8 -*-
'''
# ds1422.py      Sept 2025
# [2D] LINEAR DYNAMICAL SYSTEMS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/ds1422.pdf

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
    dx = a00*x + a01*y
    dy = a10*x + a11*y
    return [dx, dy]  

#%% INPUTS
# X and Y limits for plots / time steps 
L1 = -6; L2 = 6; nL = 999; nT = 999

# Initial values: a, x0, y0    comment / uncomment for different inputs
N = 12    # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

# a00, a01, a10, a11 = 2,1,1,2
# x0[0], y0[0] = -1, 2
# x0[1], y0[1] = -2.90, 1.31 
# x0[2], y0[2] = -1.5,-0.8
# x0[3], y0[3] = 1.00, -2.00
# x0[4], y0[4] = 3.00, -1.50
# x0[5], y0[5] = 1.50, -0.50
# x0[6], y0[6] = 0.10, 0.10
# t2 = 5
# num = 7

a00, a01, a10, a11 = -0.3,0.4,-0.2,0.3
x0[0], y0[0] = -1, 2
x0[1], y0[1] = -2.90, 1.31 
x0[2], y0[2] = -1.5,-0.8
x0[3], y0[3] = 1.00, -2.00
x0[4], y0[4] = 3.00, -1.50
x0[5], y0[5] = 1.50, -0.50
x0[6], y0[6] = 0.10, 0.10
t2 = 25
num = 7   

# a00, a01, a10, a11 = 1,1,4,-2
# x0[0], y0[0] = -4, 5 
# x0[1], y0[1] = -3,5 
# x0[2], y0[2] = -2, 5
# x0[3], y0[3] = -1,5 
# x0[4], y0[4] = -1.5,5
# x0[5], y0[5] = 0, 5
# x0[6], y0[6] = 1,-5
# x0[7], y0[7] = 2,-5 
# x0[8], y0[8] = 3,-5 
# x0[9], y0[9] = 4,-5
# t2 = 1.5 
# num = 10
 
#%% SETUP
# Linear system A matrix
A = zeros([2,2])
A[0,0] = a00; A[0,1] = a01; A[1,0] = a10; A[1,1] = a11
detA = np.linalg.det(A)

# Jacobian J matrix
J = A
JeigV, JeigF = eig(J)
JeigF[:,0] = JeigF[:,0]/JeigF[1,0]
JeigF[:,1] = JeigF[:,1]/JeigF[1,1]

# Solution ODE for x and y 
xS = zeros([nT,num]) ; yS = zeros([nT,num])
t = linspace(0,t2,nT)

for c in range(num):
  u0 = [x0[c],y0[c]]
  sol = odeint(lorenz, u0, t, tfirst=True)
  xS[:,c] = sol[:,0]     
  yS[:,c] = sol[:,1] 
  
# Phase Portrait: quiver and stream plots  
xP = linspace(L1,L2,15); yP = xP
xx,yy = np.meshgrid(xP,yP)
xxDot = a00*xx + a01*yy
yyDot = a10*xx + a11*yy

# Nullclines
xN = linspace(L1,L2,nL)
if a01 != 0:
   yNx = (-a00/a01)*xN
if a11 != 0:
   yNy = (-a10/a11)*xN 
   
# Manifolds
xM0 = JeigF[0,0]; yM0 = JeigF[1,0]
xM1 = JeigF[0,1]; yM1 = JeigF[1,1]
m0 = yM0/xM0   
m1 = yM1/xM1
XM = xP
YM0 = m0*XM; YM1 = m1*XM

    
#%% Console output
print(' ')
print('A matrix: a00 = %0.2f' %a00 + '  a01 = %0.2f' %a01
      + '  a10 = %0.2f ' %  a10 + '  a11 = %0.2f' %a11)
print('Determinant A = %0.5f' % detA)
print('Initial conditions (x0, y0)')
for c in range(num):
    print('  (%0.4f, ' %x0[c] + '%0.4f)' %y0[c])
print('Eigenvalues   %0.3f '% JeigV[0] + '  %0.3f ' % JeigV[1])
print('Eigenfunctions ', JeigF[0],'   ', JeigF[1]) 


#%% GRAPHICS
# FIGURE 1: x vs xDot   y vs yDot
plt.rcParams['font.size'] = 11
plt.rcParams["figure.figsize"] = (7,3)
fig1, ax = plt.subplots(nrows=1, ncols=2)
XP = linspace(L1,L2,88); YP = XP
XX,YY = np.meshgrid(XP,XP)
XXDot = a00*XX + a01*YY
YYDot = a10*XX * a11*YY

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
plt.rcParams["figure.figsize"] = (5,3.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('x',color= 'black',fontsize = 12)
ax.set_ylabel('y',color = 'black',fontsize = 12)
ax.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_xlim([L1,L2])
ax.set_ylim([L1,L2])
ax.set_aspect('equal')
 
ax.quiver(xx,yy,xxDot,yyDot)
ax.plot(xS,yS,'g',lw = 1)

if a01 != 0:
   ax.plot(xN,yNx,'r',lw = 2)
if a11 != 0:
   ax.plot(xN,yNy,'b',lw = 2)
if (a10 == 0) and (a11 == 0):
   ax.plot([0,0],[L1,L2],'r',lw = 2) 
if (a00 == 0) and (a01 == 0):
   ax.plot([L1,L2],[0,0],'b',lw = 2)    

ax.plot(x0,y0,'go', ms = 6)
ax.plot(0,0,'ko',ms = 6)

fig3.tight_layout()


#%% FIGURE 4: Phase Portrait  streamine
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3.5)
fig4, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
axes.xaxis.grid()
axes.yaxis.grid()
axes.set_xlim([L1,L2])
axes.set_ylim([L1,L2])
axes.set_aspect('equal')
axes.streamplot(xx,yy,xxDot,yyDot)
axes.plot(xS,yS,'g',lw = 1)

if a01 != 0:
   axes.plot(xN,yNx,'r',lw = 1)
if a11 != 0:
   axes.plot(xN,yNy,'b',lw = 1)
if (a10 == 0) and (a11 == 0):
   axes.plot([0,0],[L1,L2],'r',lw = 1) 
if (a00 == 0) and (a01 == 0):
   axes.plot([L1,L2],[0,0],'b',lw = 1)    
axes.plot(x0,y0,'go', ms = 6)
axes.plot(0,0,'ko',ms = 6)
axes.plot(x0,y0,'go', ms = 6)
axes.plot(XM,YM0,'k',lw = 2)
axes.plot(XM,YM1,'k',lw = 2)
fig4.tight_layout()


#%% FIGURE 5: Eigenfunctions and manifolds
# Manifolds
xM0 = JeigF[0,0]; yM0 = JeigF[1,0]
xM1 = JeigF[0,1]; yM1 = JeigF[1,1]
m0 = yM0/xM0   
m1 = yM1/xM1
XM = xP
YM0 = m0*XM; YM1 = m1*XM


plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3.5)
fig5, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
axes.set_xlim([L1,L2])
axes.set_ylim([L1,L2])
axes.xaxis.grid()
axes.yaxis.grid()
axes.set_aspect('equal')

axes.streamplot(xx,yy,xxDot,yyDot)
axes.plot(XM,YM0,'r',lw = 2)
axes.plot(XM,YM1,'m',lw = 2)

xA = array([-3,-3,3,3])
for c in range(4):
    yA = m0*xA[c]    
    u0 = [xA[c],yA]
    sol = odeint(lorenz, u0, t, tfirst=True)
    sx = sol[:,0]     
    sy = sol[:,1] 
    axes.plot(sx,sy,'g',lw = 3)
    axes.plot(sx[0],sy[0],'go',ms=6)
    yA = m1*xA[c]    
    u0 = [xA[c],yA]
    sol = odeint(lorenz, u0, t, tfirst=True)
    sx = sol[:,0]     
    sy = sol[:,1] 
    axes.plot(sx,sy,'g',lw = 3)
    axes.plot(sx[0],sy[0],'go',ms=6)
    
axes.plot(0,0,'ko',ms = 6)    
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