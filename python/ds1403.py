# -*- coding: utf-8 -*-
'''
ds25L1401.py
Ian Cooper
Sep 25

[2D] DYNAMICAL SYSTEMS:
   PLANAR LINEAR SYSTEMS
   MASS - SPRING SYSTEM

https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation

https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds1401.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, real, imag, sqrt, exp 
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
a00 = 2; a01 = 1; a10 = 1; a11 = 2

x0, y0 = 5, 5

t1, t2 = 0, 30

nT = 9999

#%% SETUP
u0 = [x0, y0]
t = linspace(t1,t2,nT)

# Solve ODE: Fig 1
sol = odeint(lorenz, u0, t, tfirst=True)
x = sol[:,0]     
y = sol[:,1] 

# Nullclines
xN = linspace(-6,6,999)
yNx = -xN
yNy = 2*xN 

#%%  Fig 6 Fig 7: LINEAR ALGEBRA SOLUTION
A = array([[a00,a01],[a10,a11]])
L, F = eig(A)
X = np.transpose(u0)
c = np.linalg.solve(F,X)

xA = c[0]*exp(L[0]*t)*F[0,0] + c[1]*exp(L[1]*t)*F[0,1]
yA = c[0]*exp(L[0]*t)*F[1,0] + c[1]*exp(L[1]*t)*F[1,1]

# c1 = 0
x2 = c[0]*exp(L[0]*10)*F[0,0]
y2 = c[0]*exp(L[0]*10)*F[1,0]
m2 = y2/x2
yL2 = m2*xN
# c0 = 0
x1 = c[1]*exp(L[1]*10)*F[0,1]
y1 = c[1]*exp(L[1]*10)*F[1,1]
m1 = y1/x1
yL1 = m1*xN 

#%% Fig 1: t vs x
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)

txtx = 't'; txty = 'x and y'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.grid()
xP = t; yP = x 
ax.plot(xP,yP,'b',lw = 2,label = 'x')
xP = t; yP = y 
ax.plot(xP,yP,'r',lw = 2,label = 'y')
ax.legend()
fig1.tight_layout()
fig1.savefig('a1.png')

#%% Fig 3: x vs v
fig3, ax = plt.subplots(nrows=1, ncols=1)
txtx = 'x'; txty = 'y'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.grid()
xP = x; yP = y
ax.plot(xP,yP,'b',lw = 2)
ax.plot(x0,y0,'go',ms = 6)
ax.plot(0,0,'ro',ms = 6)
fig3.tight_layout()
fig3.savefig('a3.png')

#%%  Fig 4 Fig 5: Phase Portrait: quiver and stream plots  
plt.rcParams["figure.figsize"] = (5,4)
d = 6
xG = linspace(-d,d,10); yG = linspace(-d,d,10)
xx,yy = np.meshgrid(xG,yG)
xxDot = a00*xx + a01*yy
yyDot = a10*xx + a11*yy 

fig4, ax = plt.subplots(nrows=1, ncols=1)
txtx = 'x'; txty = 'y'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.text(1.8,1,'c$_1$ = 0',bbox=dict(facecolor='yellow',
                                   edgecolor = 'white',alpha=0.5))
ax.text(-3.5,4.2,'c$_0$ = 0',bbox=dict(facecolor='yellow',
                                   edgecolor = 'white',alpha=0.5))

ax.text(-2,-5,'y$_{null}$ = 0',bbox=dict(facecolor='red',
                                   edgecolor = 'white',alpha=0.5))
ax.text(-4.5,1,'x$_{null}$ = 0',bbox=dict(facecolor='red',
                                   edgecolor = 'white',alpha=0.5))

ax.plot(0,0,'ro',ms = 8)
ax.streamplot(xx,yy,xxDot,yyDot)
ax.plot(xP,yP,'b',lw = 2)

ax.plot(xN,yNx,'r')
ax.plot(xN,yNy,'r')
ax.plot(xN,yL1,'k',lw = 2)
ax.plot(xN,yL2,'k',lw = 2)

ax.axis('square')
ax.set_xlim([-6,6]); ax.set_ylim([-6,6])
ax.set_xticks
fig4.tight_layout()
fig4.savefig('a4.png')

fig5, ax = plt.subplots(nrows=1, ncols=1)
txtx = 'x'; txty = 'v'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
fig5.tight_layout()
ax.quiver(xx,yy,xxDot,yyDot)
ax.axis('square')
ax.set_xlim([-6,6]); ax.set_ylim([-6,6])
ax.set_xticks
ax.plot(0,0,'ro',ms = 6)
ax.plot(xP,yP,'b',lw = 2)
fig5.savefig('a5.png')

#%%  Fig 6 Fig 7: LINEAR ALGEBRA SOLUTION
A = array([[a00,a01],[a10,a11]])
L, F = eig(A)
X = np.transpose(u0)
c = np.linalg.solve(F,X)

xA = c[0]*exp(L[0]*t)*F[0,0] + c[1]*exp(L[1]*t)*F[0,1]
yA = c[0]*exp(L[0]*t)*F[1,0] + c[1]*exp(L[1]*t)*F[1,1]

x1 = c[1]*exp(L[1]*10)*F[0,1]
y1 = c[1]*exp(L[1]*10)*F[1,1]
m1 = y1/x1

print(' ')
print('Initial conditions', np.round(u0,3))
print('Eigenvalues', np.round(L,2))
print('Eigenvectors', np.round(F,2))
print('Eigenvalues', np.round(L,2))
print('coeff c =', np.round(c,3))

fig6, ax = plt.subplots(nrows=1, ncols=1)
txtx = 't'; txty = 'x'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.grid()
xP = t; yP = xA 
ax.plot(xP, yP,'b',lw = 2)
xP = t; yP = yA 
ax.plot(xP, yP,'r',lw = 2)
fig6.tight_layout()
fig6.savefig('a6.png')



#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
