# -*- coding: utf-8 -*-

'''
ds1510.py      OCt 2025
NONLINEAR [2D] DYNAMICAL SYSTEMS
HOMOCLINIC BIFURCATIONS

Ian Cooper: matlabvisualphysics@gmail.com

   https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
   https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS1511.pdf

'''


# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')
tStart = time.time()

# FUNCTIONS  Solve ODE for x,y    x = R   y = theta 
def lorenz(t, state):    
    x, y = state
    dx = y
    dy = r*y + x - x**2 + x*y
    return [dx, dy]  

#%%
r  = 1
x0 = 0.9
y0 = -0.50
# Solution ODE for x and y 
t1 = 0; t2 = 1; nT = 999
t = linspace(t1,t2,nT)
u0 = [x0,y0]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]            


yC = (x0**2 - x0**6/3)**0.5

#%%
rJ = r
# Jacobian matrix and eigenvalues
J0 = np.array([[0,1],[1,rJ]])
ev0, J0 = eig(J0)
print(ev0)
JP = np.array([[0,1],[-1,rJ+1]])
evP, efP1 = eig(JP)
print(evP)
JM = np.array([[0,1],[3,rJ-1]])
evM, efPM = eig(JM)
print(evM)



#%%  GRAPHICS  Phase portrait orbits
# Figure 1
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3.4,3.2)
fig1, axes = plt.subplots(nrows=1, ncols=1)
#fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
#                    right = 0.95, hspace = 0.5,wspace=0.40)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
#axes.set_title('x(0) = %2.3f' % x0 + '  y(0) = %2.3f' % y0, fontsize = 10)
axes.set_title('r = %2.3f' % r, fontsize = 12)

# axes.set_xlim([-1.5, 1.7])
# axes.set_xticks(np.arange(-1.5,1.6,0.5))
#axes.set_ylim([-1.5, 1.5])
# axes.set_yticks(np.arange(-1,1.1,0.5))

axes.xaxis.grid()
axes.yaxis.grid()

axes.plot(xS[0], yS[0],'go',ms = 8)
axes.plot(xS, yS,'b',lw = 2)
axes.plot(0, 0,'ko',ms = 6)
axes.plot(1, 0,'ko',ms = 6)
axes.plot(-1, 0,'ko',ms = 6)

#axes.set_aspect('equal', 'box')
fig1.tight_layout()

fig1.savefig('a1.png')

#%%
# #%% Phase Portrait quiver plot  
# # Figure 2
# x1 = -1; x2 = 1; nX = 21
# xP = linspace(x1,x2,nX)
# y1 = -1; y2 = 1; nY = 21
# yP = linspace(y1,y2,nX)
# xx,yy = np.meshgrid(xP,yP)
# xxDot = r*xx + xx**3 - xx**5
# yyDot = -yy

# plt.rcParams['font.size'] = 10
# plt.rcParams["figure.figsize"] = (4,4)
# fig1, axes = plt.subplots(nrows=1, ncols=1)
# fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
#                      right = 0.95, hspace = 0.5,wspace=0.40)
  
# axes.set_xlabel('x',color= 'black',fontsize = 12)
# axes.set_ylabel('y',color = 'black',fontsize = 12)
# axes.set_title('r = %2.1f' % r, fontsize = 14)
# axes.set_xlim([-1.1, 1.1])
# axes.set_yticks(np.arange(-1,1.1,0.5))
# axes.set_ylim([-1.1, 1.1])
# axes.set_yticks(np.arange(-1,1.1,0.5))
# axes.plot(0,0,'ro',ms = 6)
# axes.streamplot(xx,yy,xxDot,yyDot, density = 0.8, color = 'blue')

# fig1.tight_layout()
# axes.set_aspect('equal', 'box')

# fig1.savefig('a2.png')


# #%%  TIME EVOLUTION PLOTS  T r x y
# Figure 3
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3.6)
fig1, axes = plt.subplots(nrows=1, ncols=2)
# fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.15,\
#                     right = 0.95, hspace = 0.5,wspace=0.5)
#fig1.suptitle('r = %2.1f' % r, fontsize = 14)  
  
Cg = 0   
#axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
#axes[Rg,Cg].set_ylabel(r'$\theta$',color = 'black',fontsize = 12)
# #axes[C].set_xlim([-5, 5])
# #axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[Rg,Cg].set_ylim([0, 26])
# #axes[R,C].set_yticks(np.arange(-20,81,20))
axes[Cg].xaxis.grid()
axes[Cg].yaxis.grid()
axes[Cg].plot(t, xS,'b',lw = 2)

Rg = 0;Cg = 1   
axes[Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Cg].set_ylabel('R',color = 'black',fontsize = 12)
# #axes[C].set_xlim([-5, 5])
# #axes[C].set_xticks(np.arange(-5,5.1,1))
# #axes[R,C].set_ylim([-200, 200])
# #axes[R,C].set_yticks(np.arange(-20,81,20))
axes[Cg].xaxis.grid()
axes[Cg].yaxis.grid()
axes[Cg].plot(t, yS,'b',lw = 2)

# Rg = 1;Cg = 0   
# axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
# axes[Rg,Cg].set_ylabel('x',color = 'black',fontsize = 12)
# # #axes[C].set_xlim([-5, 5])
# # #axes[C].set_xticks(np.arange(-5,5.1,1))
# # #axes[R,C].set_ylim([-200, 200])
# # #axes[R,C].set_yticks(np.arange(-20,81,20))
# axes[Rg,Cg].xaxis.grid()
# axes[Rg,Cg].yaxis.grid()
# axes[Rg,Cg].plot(t, xS,'b',lw = 2)

# Rg = 1;Cg = 1   
# axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
# axes[Rg,Cg].set_ylabel('y',color = 'black',fontsize = 12)
# # #axes[C].set_xlim([-5, 5])
# # #axes[C].set_xticks(np.arange(-5,5.1,1))
# # #axes[R,C].set_ylim([-200, 200])
# # #axes[R,C].set_yticks(np.arange(-20,81,20))
# axes[Rg,Cg].xaxis.grid()
# axes[Rg,Cg].yaxis.grid()
# axes[Rg,Cg].plot(t, yS,'b',lw = 2)

fig1.tight_layout()

# fig1.savefig('a3.png')


# #%%  Jacobian matrix and eigenvalues
# J = np.array([[r,-1],[1,r]])
# Jev, Jef = eig(J)
# print(J)
# print(Jev)


# #%%   Rdot as a function of r
# # Figure 4
# R1 = 0; R2 = 1.2; N = 599
# rp = r
# Rp = linspace(R1,R2,N)
# Rdot = rp*Rp + Rp**3 - Rp**5

# RP = ((1+(1+4*rp)**0.5)/2)**0.5

# fRP = r + 3*RP**2 -5*RP**4

# print(RP,fRP)


# plt.rcParams['font.size'] = 10
# plt.rcParams["figure.figsize"] = (4,3)
# fig1, axes = plt.subplots(nrows=1, ncols=1)
# #fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
# #                    right = 0.95, hspace = 0.5,wspace=0.40)
  
# axes.set_xlabel('R',color= 'black',fontsize = 12)
# axes.set_ylabel('Rdot',color = 'black',fontsize = 12)
# axes.set_title('r = %2.1f' % rp, fontsize = 14)
# # axes.set_xlim([-1.1, 1.1])
# # axes.set_yticks(np.arange(-1,1.1,0.5))
# # axes.set_ylim([-1.1, 1.1])
# # axes.set_yticks(np.arange(-1,1.1,0.5))
# axes.xaxis.grid()
# axes.yaxis.grid()

# axes.plot(0,0,'ro',ms = 8)
# axes.plot(RP,0,'bo',ms = 8)
# axes.plot(Rp,Rdot,'b',lw = 2)

# fig1.tight_layout()

# fig1.savefig('a4.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
