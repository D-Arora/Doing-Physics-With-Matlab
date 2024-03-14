# -*- coding: utf-8 -*-

# cs120.py    Feb 2024

# NONLINEAR [2D] DYNAMICAL SYSTEMS
# FIXED POINTS, STABILITY ANALYSIS, BIFURCATIONS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs120.pdf


# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

tStart = time.time()

#%% FUNCTIONS  Solve ODE for x,y    x = theta   y = omega 
def lorenz(t, state):    
    x, y = state
    dx = r - x**2
    dy = -y
    return [dx, dy]  

r = -9

# x vs xDot   y vs yDot -------------------------------------------------------
y1 = -10; y2 = 10; nY = 599
y = linspace(y1,y2,nY)
yDot = - y

x1 = -5; x2 = 5; nX = 599
x = linspace(x1,x2,nX)
xDot = r - x**2


# Solution ODE for x and y  ---------------------------------------------------
t1 = 0; t2 = 0.6; nT = 999
t = linspace(t1,t2,nT)
u0 = [0.8,0.01]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]            


# Phase Portrait quiver plot  -------------------------------------------------
x1 = -5; x2 = 5; nX = 9
xP = linspace(x1,x2,nX)
y1 = -5; y2 = 5; nY = 9
yP = linspace(y1,y2,nX)

xx,yy = np.meshgrid(xP,yP)
xxDot = r - xx**2
yyDot = -yy

# Jacobian matrix and eigenvalues
J = np.array([[0,0],[0,-1]])
Jev, Jef = eig(J)
print(Jev)


#%% GRAPHICS
# FIGURE 1: x vs xDot   y vs yDot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,3)
fig1, axes = plt.subplots(nrows=1, ncols=2)
fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
                    right = 0.95, hspace = 0.5,wspace=0.40)

    
C = 0   
axes[C].set_xlabel('x',color= 'black',fontsize = 12)
axes[C].set_ylabel('$x_{dot}$',color = 'black',fontsize = 12)
axes[C].set_title('r = %2.1f' % r, fontsize = 14)
axes[C].set_xlim([-5, 5])
axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(x, xDot,'b',lw = 2)
#axes[C].plot(-np.sqrt(r),0,'ro',ms = 6)
#axes[C].plot(np.sqrt(r),0,'bo',ms = 6)

C = 1   
axes[C].set_xlabel('y',color= 'black',fontsize = 12)
axes[C].set_ylabel('$y_{dot}$',color = 'black',fontsize = 12)
axes[C].set_title('r = %2.1f' % r, fontsize = 14)
#axes[R,C].set_xlim([-6, 6])
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_xticks(np.arange(-6,7,2))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(y, yDot,'b',lw = 2)
axes[C].plot(0,0,'bo',ms = 6)

fig1.savefig('a1.png')

# FIGURE 2: t vs x   t vs y
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,3)
fig1, axes = plt.subplots(nrows=1, ncols=2)
fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.15,\
                    right = 0.95, hspace = 0.5,wspace=0.40)

    
C = 0   
axes[C].set_xlabel('t',color= 'black',fontsize = 12)
axes[C].set_ylabel('x',color = 'black',fontsize = 12)
axes[C].set_title('r = %2.1f' % r, fontsize = 14)
#axes[C].set_xlim([-5, 5])
#axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(t, xS,'b',lw = 2)

C = 1   
axes[C].set_xlabel('t',color= 'black',fontsize = 12)
axes[C].set_ylabel('y',color = 'black',fontsize = 12)
axes[C].set_title('r = %2.1f' % r, fontsize = 14)
#axes[R,C].set_xlim([-6, 6])
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_xticks(np.arange(-6,7,2))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[C].xaxis.grid()
axes[C].yaxis.grid()
axes[C].plot(t, yS,'b',lw = 2)
fig1.savefig('a2.png')


# FIGURE 3: x vs y   Phase Portrait   /  nullclines
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3,3)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.18,\
                    right = 0.95, hspace = 0.5,wspace=0.40)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)
#axes.set_xlim([-20, 0])
#axes.set_xticks(np.arange(-20,0.1,5))
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes.xaxis.grid()
axes.yaxis.grid()
axes.plot(xS, yS,'b',lw = 2)
axes.plot(xS[0], yS[0],'go',ms = 8)
#axes.plot([-20,0], [0,0],'r',lw = 2)
#axes.plot([0,0.0], [0,0],'r',lw = 2)
axes.plot([0,0.8], [0,0],'r',lw = 2)
axes.plot([0,0], [0,0.010],'m',lw = 3)
fig1.savefig('a3.png')


# FIGURE 4: Phase Portrait  quiver plot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
                    right = 0.95, hspace = 0.5,wspace=0.40)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)
#axes[C].set_xlim([-5, 5])
#axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes.xaxis.grid()
axes.yaxis.grid()
#axes.quiver(xx,yy,xxDot,yyDot)
axes.streamplot(xx,yy,xxDot,yyDot)
fig1.savefig('a4.png')

# FIGURE 5: Phase Portrait  quiver plot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
                    right = 0.95, hspace = 0.5,wspace=0.40)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)
#axes[C].set_xlim([-5, 5])
#axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes.xaxis.grid()
axes.yaxis.grid()
axes.quiver(xx,yy,xxDot,yyDot)
#axes.streamplot(xx,yy,xxDot,yyDot)
fig1.savefig('a5.png')


# FIGURE 6: 
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,3)
fig1, axes = plt.subplots(nrows=1, ncols=2)
fig1.subplots_adjust(top = 0.90, bottom = 0.18, left = 0.13,\
                    right = 0.95, hspace = 0.5,wspace=0.40)
fig1.suptitle('r = %2.1f' % r, fontsize = 14)
    
C = 0   
axes[0].set_xlabel('x',color= 'black',fontsize = 12)
axes[0].set_ylabel('y',color = 'black',fontsize = 12)
#axes[0].set_title('r = %2.1f' % r, fontsize = 14)
#axes[C].set_xlim([-5, 5])
#axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[0].xaxis.grid()
axes[0].yaxis.grid()
axes[0].quiver(xx,yy,xxDot,yyDot)
#axes[0].plot(-3,0,'ro',ms = 8)
#axes[0].plot(3,0,'bo',ms = 8)

C = 1   
axes[1].set_xlabel('x',color= 'black',fontsize = 12)
axes[1].set_ylabel('y',color = 'black',fontsize = 12)
#axes[1].set_title('r = %2.1f' % r, fontsize = 14)
#axes[C].set_xlim([-5, 5])
#axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[R,C].set_ylim([-200, 200])
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[1].xaxis.grid()
axes[1].yaxis.grid()
#axes[1].plot(-3,0,'ro',ms = 8)
#axes[1].plot(3,0,'bo',ms = 8)
axes[1].streamplot(xx,yy,xxDot,yyDot)

fig1.savefig('a6.png')

tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
