# -*- coding: utf-8 -*-
'''
ds120.py      Oct 2025
NONLINEAR [2D] DYNAMICAL SYSTEMS
SUPERCRITICAL HOPF BIFURCATION

Ian Cooper: matlabvisualphysics@gmail.com

   https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
   https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS1509.pdf

'''

# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, array, zeros 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')
tStart = time.time()

# FUNCTIONS  Solve ODE for x,y    x = R   y = theta 
def lorenz(t, state):    
    x, y = state
    dx = r*x - x**3 
    dy = 2*pi
    return [dx, dy]  

#%%  Jacobian matrix and eigenvalues
r = 1

J = np.array([[r,1],[-1,r]])
Jev, Jef = eig(J)
print(J)
print(Jev)

#%% Figure 1: r vs Rdot
R = linspace(0,2.2,599)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4.5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.grid()
ax.set_xlabel('R',color= 'black',fontsize = 12)
ax.set_ylabel('R$_{dot}$')
r = -4; Rdot = r*R - R**3; ax.plot(R,Rdot,'r',lw = 2,label = 'r = -4')
r = 0;  Rdot = r*R - R**3; ax.plot(R,Rdot,'b',lw = 2,label = 'r = 0')
r = 4;  Rdot = r*R - R**3; ax.plot(R,Rdot,'m',lw = 2,label = 'r = +4')
ax.legend()
fig1.tight_layout()
fig1.savefig('a1.png')


#%%  Solution ODE for x and y
R0, T0 = 4.8,1.5*pi/4
r = 4
t1 = 0; t2 = 10; nT = 9999
t2 = 10
t = linspace(t1,t2,nT)
u0 = [R0,T0]
sol = odeint(lorenz, u0, t, tfirst=True)
R = sol[:,0]     
T = sol[:,1]            
xS = R*cos(T)
yS = R*sin(T)


#%% FIGURE 2: Time evoltion plots  t vs R
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4.5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x, y')
ax.set_title('r = %2.3f' % r, fontsize = 14)
ax.grid()
ax.plot(t, xS,'b',lw = 2)
ax.plot(t, yS,'r',lw = 2)
fig2.tight_layout()
fig2.savefig('a2.png')


#%% FIGURE 3: Phase portrait x vs y
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x, y')
ax.set_title('r = %2.2f' % r, fontsize = 14)

if r <= 0:
   ax.set_xlim([-1.1, 1.1])
   ax.set_xticks(np.arange(-1,1.1,0.5))
   ax.set_ylim([-1.1, 1.1])
   ax.set_yticks(np.arange(-1,1.1,0.5))
ax.grid()

ax.plot(xS, yS,'b',lw = 2)
ax.plot(xS[0], yS[0],'go',lw = 2)

ax.set_box_aspect(1)
fig3.tight_layout()
fig3.savefig('a3.png')


#%% FIGURE 4: Phase Portrait streamplot  
x1 = -5; x2 = 5; nX = 21
xP = linspace(x1,x2,nX)
y1 = -5; y2 = 5; nY = 21
yP = linspace(y1,y2,nX)
xx,yy = np.meshgrid(xP,yP)
xxDot = r*xx - xx**3
yyDot = -yy

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,3)
fig4, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.1f' % r, fontsize = 14)

axes.plot(0,0,'ro',ms = 6)
axes.streamplot(xx,yy,xxDot,yyDot, linewidth = 1, density = 0.8,
    color = 'b')

axes.plot(xS, yS,'g',lw = 2)

fig4.tight_layout()
axes.set_aspect('equal', 'box')
fig4.tight_layout()
fig4.savefig('a4.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
