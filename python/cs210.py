# -*- coding: utf-8 -*-
'''
# cs200.py      June 2025
# NONLINEAR [2D] DYNAMICAL SYSTEMS
# FIXED POINTS, STABILITY ANALYSIS, BIFURCATIONS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs200.pdf

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
    dx = x
    dy = x**2 + y**2 - 1
    return [dx, dy]  

#%% INPUTS
# X and Y limits for plots 
L1 = -3; L2 = 3; nL = 599

# Time span 
t1 = 0; t2 = 1; nT = 599

# Initial vales x0, y0    comment / uncomment for different inputs
N = 11    # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

# manual IC
x0[0] = 0.2; y0[0] = 1.2
x0[1] = 0; y0[1] = 1.2
x0[2] = -0.2; y0[2] = 1.2
x0[3] = 0.5; y0[3] = 0.8
x0[4] = -0.5; y0[4] = 0.8
x0[5] = 0.5; y0[5] = 0.4
x0[6] = -0.5; y0[6] = 0.4
x0[7] = 0.5; y0[7] = -0.2
x0[8] = -0.5; y0[8] = -0.2
x0[9] = 0.5; y0[9] = -1.2
x0[10] = -0.5; y0[10] = -1.2


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
xP = linspace(L1,L2,15); yP = xP
xx,yy = np.meshgrid(xP,yP)
xxDot = xx
yyDot = xx**2 + yy**2 - 1


# Nullclines X [0] / Y [1]
xN0 = zeros(nL); yN0 = linspace(L1,L2,nL) 
xN1 = linspace(-1,1,nL)
yN1 = sqrt(1-xN1**2)

# Critical points
xC = array([0, 0]); yC = array([1,-1])

# Jacobian matrix J
J11 = array([1,1]); J12 = ([0,0]); J21 = 2*xC; J22 = 2*yC

J0 = array([[J11[0],J12[0]],[J21[0],J22[0]]])
J1 = array([[J11[1],J12[1]],[J21[1],J22[1]]])

# Eignevalues V, eigenvector (eigenfunction) F
V0, F0 = eig(J0); V1,F1 = eig(J1)

    
#%% Console output
print(' ')
print('Critical point 0(xC, yC)')
print('   (%0.0f  ' % xC[0] + ', %0.3f)' %yC[0] ) 
print('Jacobian matrix J0')
print(J0)
print('Eigenvalues J0' )
print('   (%0.3f  ' % V0[0] + ', %0.3f)' %V0[1] ) 
print('Eigenfunctions J0')
print(F0)

print('Critical point 1(xC, yC)')
print('   (%0.0f  ' % xC[1] + ', %0.3f)' %yC[1] ) 
print('Jacobian matrix J1')
print(J1)
print('Eigenvalues J1' )
print('   (%0.3f  ' % V1[0] + ', %0.3f)' %V1[1] ) 
print('Eigenfunctions J1')
print(F1)



#%% GRAPHICS
# FIGURE 1: x vs xDot   y vs yDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,4)
fig1, ax = plt.subplots(nrows=1, ncols=2)
XP = linspace(L1,L2,88); YP = XP
XX,YY = np.meshgrid(XP,XP)
XXDot = XX
YYDot = XX**2 + YY**2 - 1

C = 0   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('x$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([L1,L2])
cf = ax[C].pcolor(XX,YY,XXDot, cmap='jet') 
ax[C].plot([0,0],[L1,L2],'w', lw = 2)
ax[C].set_box_aspect(1)
fig1.colorbar(cf, ax=ax[C])

C = 1   
ax[C].set_xlabel('x',color= 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].set_title('y$_{dot}$',color = 'black',fontsize = 14)
ax[C].set_xlim([L1,L2])
cf = ax[C].pcolor(XX,YY,YYDot, cmap='jet') 
fig1.colorbar(cf, ax=ax[C])
ax[C].set_xticks([-2,-1,0,1,2])
ax[C].set_yticks([-2,-1,0,1,2])
ax[C].plot(xN1,yN1,'w', lw = 2)
ax[C].plot(xN1,-yN1,'w', lw = 2)
ax[C].set_box_aspect(1)

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

axes.plot(xC,yC,'ro', ms = 7)

axes.plot(xN0,yN0,'r', lw = 1)
axes.plot(xN1,yN1,'b', lw = 1)
axes.plot(xN1,-yN1,'b', lw = 1)
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



axes.plot(xC,yC,'ro', ms = 7)

axes.plot(xN0,yN0,'r', lw = 1)
axes.plot(xN1,yN1,'b', lw = 1)
axes.plot(xN1,-yN1,'b', lw = 1)

axes.plot(xS,yS,'g',lw = 3)
axes.plot(x0,y0,'go', ms = 6)

axes.set_box_aspect(1)

fig5.tight_layout()

#%%   FIGURE 3: slope plot dy/dx
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)
fig3, ax = plt.subplots(nrows=1, ncols=1)
XP = linspace(-2,2,255); YP = XP
XX,YY = np.meshgrid(XP,XP)
ZZ = (XX**2+YY**2-1)/(XX+1e-16)
ZZ[(ZZ)>20]=20;ZZ[(ZZ)<-20] = -20  
ax.set_xlabel('x',color= 'black',fontsize = 12)
ax.set_ylabel('y',color = 'black',fontsize = 12)
ax.set_title('slope function  dy/dx',color = 'black',fontsize = 14)
ax.set_xticks([-2,-1,0,1,2])
ax.set_yticks([-2,-1,0,1,2])
cf = ax.pcolor(XX,YY,abs(ZZ)**0.3, cmap='jet') 
ax.set_box_aspect(1)
fig3.colorbar(cf, ax=ax)
fig3.tight_layout()


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