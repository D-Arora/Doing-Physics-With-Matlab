# -*- coding: utf-8 -*-

# cs130.py                 Mar 2024

# NONLINEAR [2D] DYNAMICAL SYSTEMS
# FIXED POINTS, STABILITY ANALYSIS, BIFURCATIONS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs120.pdf

# Subcritical Hopf Bifurcation



# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

tStart = time.time()

# FUNCTIONS  Solve ODE for x,y    x = R   y = theta 
def lorenz(t, state):    
    x, y = state
    dx = r*x + x**3 - x**5
    dy = 1
    return [dx, dy]  


#%%
# Bifurcation parameter
r = -0.5

# Input initial conditions  radius R and polar angle T (theta) >>>>>
T0 = 0*pi
#R0 = np.sqrt((1-np.sqrt(1+4*r))/2)
#print(R0)
R0 = 1

# Solution ODE for x and y 
t1 = 0; t2 = 15; nT = 999
t = linspace(t1,t2,nT)
u0 = [R0,T0]
sol = odeint(lorenz, u0, t, tfirst=True)
R = sol[:,0]     
T = sol[:,1]            
xS = R*cos(T)
yS = R*sin(T)

#%%  GRAPHICS  Phase portrait orbits
# Figure 1
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3.6,3)
fig1, axes = plt.subplots(nrows=1, ncols=1)
#fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
#                    right = 0.95, hspace = 0.5,wspace=0.40)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)
#axes.set_xlim([-0.1, 0.1])
# axes.set_yticks(np.arange(-1,1.1,0.5))
# axes.set_ylim([-1.1, 1.1])
# axes.set_yticks(np.arange(-1,1.1,0.5))
axes.xaxis.grid()
axes.yaxis.grid()

axes.plot(xS[0], yS[0],'go',ms = 8)
axes.plot(xS, yS,'b',lw = 2)

#axes.set_aspect('equal', 'box')
fig1.tight_layout()

fig1.savefig('a1.png')


#%% Phase Portrait quiver plot  
# Figure 2
x1 = -1; x2 = 1; nX = 21
xP = linspace(x1,x2,nX)
y1 = -1; y2 = 1; nY = 21
yP = linspace(y1,y2,nX)
xx,yy = np.meshgrid(xP,yP)
xxDot = r*xx + xx**3 - xx**5
yyDot = -yy

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3.6,3.2)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
                     right = 0.95, hspace = 0.5,wspace=0.40)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)
axes.set_xlim([-1.1, 1.1])
axes.set_yticks(np.arange(-1,1.1,0.5))
axes.set_ylim([-1.1, 1.1])
axes.set_yticks(np.arange(-1,1.1,0.5))
axes.plot(0,0,'bo',ms = 6)
axes.streamplot(xx,yy,xxDot,yyDot, density = 0.8, color = 'blue')

fig1.tight_layout()
axes.set_aspect('equal', 'box')

fig1.savefig('a2.png')


#%%  TIME EVOLUTION PLOTS  T r x y
# Figure 3
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6,6)
fig1, axes = plt.subplots(nrows=2, ncols=2)
# fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.15,\
#                     right = 0.95, hspace = 0.5,wspace=0.5)
fig1.suptitle('r = %2.1f' % r, fontsize = 14)  
  
Rg = 0;Cg = 0   
axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Rg,Cg].set_ylabel(r'$\theta$',color = 'black',fontsize = 12)
# #axes[C].set_xlim([-5, 5])
# #axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[Rg,Cg].set_ylim([0, 26])
# #axes[R,C].set_yticks(np.arange(-20,81,20))
axes[Rg,Cg].xaxis.grid()
axes[Rg,Cg].yaxis.grid()
axes[Rg,Cg].plot(t, T,'b',lw = 2)

Rg = 0;Cg = 1   
axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Rg,Cg].set_ylabel('R',color = 'black',fontsize = 12)
# #axes[C].set_xlim([-5, 5])
# #axes[C].set_xticks(np.arange(-5,5.1,1))
# #axes[R,C].set_ylim([-200, 200])
# #axes[R,C].set_yticks(np.arange(-20,81,20))
axes[Rg,Cg].xaxis.grid()
axes[Rg,Cg].yaxis.grid()
axes[Rg,Cg].plot(t, R,'b',lw = 2)

Rg = 1;Cg = 0   
axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Rg,Cg].set_ylabel('x',color = 'black',fontsize = 12)
# #axes[C].set_xlim([-5, 5])
# #axes[C].set_xticks(np.arange(-5,5.1,1))
# #axes[R,C].set_ylim([-200, 200])
# #axes[R,C].set_yticks(np.arange(-20,81,20))
axes[Rg,Cg].xaxis.grid()
axes[Rg,Cg].yaxis.grid()
axes[Rg,Cg].plot(t, xS,'b',lw = 2)

Rg = 1;Cg = 1   
axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Rg,Cg].set_ylabel('y',color = 'black',fontsize = 12)
# #axes[C].set_xlim([-5, 5])
# #axes[C].set_xticks(np.arange(-5,5.1,1))
# #axes[R,C].set_ylim([-200, 200])
# #axes[R,C].set_yticks(np.arange(-20,81,20))
axes[Rg,Cg].xaxis.grid()
axes[Rg,Cg].yaxis.grid()
axes[Rg,Cg].plot(t, yS,'b',lw = 2)

fig1.tight_layout()

fig1.savefig('a3.png')


#%%  Jacobian matrix and eigenvalues
J = np.array([[r,-1],[1,r]])
Jev, Jef = eig(J)
print(J)
print(Jev)


#%%   Rdot as a function of r
# Figure 4
R1 = 0; R2 = 1; N = 599
rp = r
Rp = linspace(R1,R2,N)
Rdot = rp*Rp + Rp**3 - Rp**5

RP = ((1+(1+4*rp)**0.5)/2)**0.5
RM = ((1-(1+4*rp)**0.5)/2)**0.5

fRP = r + 3*RP**2 -5*RP**4
fRM = r + 3*RM**2 -5*RM**4


print(RP,fRP)
print(RM,fRM)

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,3)
fig1, axes = plt.subplots(nrows=1, ncols=1)
#fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
#                    right = 0.95, hspace = 0.5,wspace=0.40)
  
axes.set_xlabel('R',color= 'black',fontsize = 12)
axes.set_ylabel('Rdot',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % rp, fontsize = 14)
# axes.set_xlim([-1.1, 1.1])
# axes.set_yticks(np.arange(-1,1.1,0.5))
# axes.set_ylim([-1.1, 1.1])
# axes.set_yticks(np.arange(-1,1.1,0.5))
axes.xaxis.grid()
axes.yaxis.grid()

axes.plot(0,0,'bo',ms = 8)
axes.plot(RM,0,'ro',ms = 8)
axes.plot(RP,0,'bo',ms = 8)
axes.plot(Rp,Rdot,'b',lw = 2)

fig1.tight_layout()

fig1.savefig('a4.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
