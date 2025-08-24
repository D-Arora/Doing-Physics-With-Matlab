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
    dx = x*(x-1)*(2-x)*(3-x)
    dy = -1
    return [dx, dy]  

#%% INPUTS
# X and Y limits for plots 
L1 = -4; L2 = 4; nL = 599

# Time span 
t1 = 0; t2 = 2*pi; nT = 599

# Initial vales x0, y0    comment / uncomment for different inputs
N = 3    # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

# manual IC
x0[0], y0[0] = 1,0
x0[1], y0[1] = 2,pi
x0[2], y0[2] = 0.5,pi/2


#%%
x = linspace(0,4,299)
fn = x*(x-1)*(2-x)*(3-x)
plt.plot(x,fn)
# x0[3], y0[3] = -0.5,-q

# x0[4], y0[4] = 0,-q
# x0[5], y0[5] = 0.5,q
# x0[6], y0[6] = 1,-q
# x0[7], y0[7] = -2,q
# x0[8], y0[8] = -1,q
# x0[9], y0[9] = 2,q
# x0[9], y0[9] = -0.5,0
# x0[10], y0[10] = 0.5,0


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
xxDot = xx*(xx-1)*(2-xx)*(3-xx)
yyDot = -np.ones([15,15])



# # Nullclines X [0] / Y [1]
# yN0 = zeros(nL); xN0 = linspace(L1,L2,nL) 
# xN1 = linspace(L1,L2,nL)
# yN1 = xN1**3-xN1

# # Critical points
# xC = array([0, 1,-1]); yC = array([0,0,0])

# # Jacobian matrix J
# J11 = array([0,0,0]); J12 = ([1,1,1]); J21 = 1-3*xC**2; J22  = ([1,1,1])

# J0 = array([[J11[0],J12[0]],[J21[0],J22[0]]])
# J1 = array([[J11[1],J12[1]],[J21[1],J22[1]]])
# J2 = array([[J11[2],J12[2]],[J21[2],J22[2]]])
# # Eignevalues V, eigenvector (eigenfunction) F
# V0, F0 = eig(J0); V1,F1 = eig(J1); V2,F2 = eig(J2)

    


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
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (7,4)
# fig1, ax = plt.subplots(nrows=1, ncols=2)
# XP = linspace(L1,L2,88); YP = XP
# XX,YY = np.meshgrid(XP,XP)
# XXDot = YY
# YYDot = XX*(1-XX**2) + YY

# C = 0   
# ax[C].set_xlabel('x',color= 'black',fontsize = 12)
# ax[C].set_ylabel('y',color = 'black',fontsize = 12)
# ax[C].set_title('x$_{dot}$',color = 'black',fontsize = 14)
# ax[C].set_xlim([L1,L2])
# cf = ax[C].pcolor(XX,YY,XXDot, cmap='jet') 
# ax[C].plot([L1,L2],[0,0],'w', lw = 2)
# ax[C].set_box_aspect(1)
# fig1.colorbar(cf, ax=ax[C])

# C = 1   
# ax[C].set_xlabel('x',color= 'black',fontsize = 12)
# ax[C].set_ylabel('y',color = 'black',fontsize = 12)
# ax[C].set_title('y$_{dot}$',color = 'black',fontsize = 14)

# cf = ax[C].pcolor(XX,YY,YYDot, cmap='jet') 
# fig1.colorbar(cf, ax=ax[C])


# ax[C].plot(xN1,yN1,'w', lw = 2)

# ax[C].set_box_aspect(1)

# ax[C].set_xlim([L1,L2])
# ax[C].set_ylim([L1,L2])
# ax[C].set_xticks([-2,-1,0,1,2])
# ax[C].set_yticks([-2,-1,0,1,2])

# fig1.tight_layout()


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
axes[C].plot(t, yS/(2*pi),lw = 2)
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

# axes.set_box_aspect(1)

fig4.tight_layout()

# #%% FIGURE 5: Phase Portrait  streamine
# plt.rcParams['font.size'] = 10
# plt.rcParams["figure.figsize"] = (5,4)
# fig5, axes = plt.subplots(nrows=1, ncols=1)
  
# axes.set_xlabel('x',color= 'black',fontsize = 12)
# axes.set_ylabel('y',color = 'black',fontsize = 12)
# axes.set_title('$\Delta$t = %0.1f' %t2,color = 'black',fontsize = 14)
# axes.xaxis.grid()
# axes.yaxis.grid()
# axes.set_xlim([L1,L2])
# axes.set_ylim([L1,L2])
# axes.streamplot(xx,yy,xxDot,yyDot)



# axes.plot(xC,yC,'ro', ms = 7)

# axes.plot(xN0,yN0,'r', lw = 1)
# axes.plot(xN1,yN1,'b', lw = 1)

# axes.plot(xS,yS,'g',lw = 3)
# axes.plot(x0,y0,'go', ms = 6)

# axes.set_box_aspect(1)

# fig5.tight_layout()

# #%%   FIGURE 3: slope plot dy/dx
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (5,4)
# fig3, ax = plt.subplots(nrows=1, ncols=1)
# XP = linspace(2,-2,555); YP = XP
# XX,YY = np.meshgrid(XP,XP)
# ZZ = XX*(1-XX**2)/(YY+1e-16)
# ZZ[(ZZ)>20]=20;ZZ[(ZZ)<-20] = -20
  
# ax.set_xlabel('x',color= 'black',fontsize = 12)
# ax.set_ylabel('y',color = 'black',fontsize = 12)
# ax.set_title('slope function  dy/dx',color = 'black',fontsize = 14)
# ax.set_xticks([-2,-1,0,1,2])
# ax.set_yticks([-2,-1,0,1,2])
# cf = ax.pcolor(XX,YY,abs(ZZ)**0.3, cmap='jet') 
# ax.plot(xC,yC,'ko', ms = 4)
# ax.set_box_aspect(1)
# fig3.colorbar(cf, ax=ax)
# fig3.tight_layout()


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