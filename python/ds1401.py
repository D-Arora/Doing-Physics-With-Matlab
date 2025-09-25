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

#%% FUNCTIONS  Solve ODE for x,(y = v)    
def lorenz(t, state): 
    x, y = state
    dx = y
    dy = K1*x + K2*y
    return [dx, dy]  

#%% INPUTS: mass m, spring constant k and damping constant q
m, k, q = 1, 4*pi**2, 1.20

# Initial conditions / time span
x0, y0 = 0.8, 1.2

t1, t2 = 0, 5

nT = 9999

#%% SETUP
K1 = -k/m
K2 = -q/m
w = sqrt(-K1)
T = 2*pi/w
u0 = [x0, y0]
t = linspace(t1,t2,nT)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)

# Solve ODE: Fig 1
sol = odeint(lorenz, u0, t, tfirst=True)
x = sol[:,0]     
v = sol[:,1] 

#%% Fig 1: t vs x
fig1, ax = plt.subplots(nrows=1, ncols=1)

txtx = 't'; txty = 'x'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.grid()
xP = t; yP = x 
ax.plot(xP,yP,'b',lw = 2)
fig1.tight_layout()

#%% Fig 2: t vs v
fig2, ax = plt.subplots(nrows=1, ncols=1)
txtx = 't'; txty = 'v'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.grid()
xP = t; yP = v
ax.plot(xP,yP,'b',lw = 2)
fig2.tight_layout()

#%% Fig 3: x vs v
fig3, ax = plt.subplots(nrows=1, ncols=1)
txtx = 'x'; txty = 'v'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.grid()
xP = x; yP = v
ax.plot(xP,yP,'b',lw = 2)
ax.plot(x0,y0,'go',ms = 6)
ax.plot(0,0,'ro',ms = 6)
fig3.tight_layout()

#%%  Fig 4 Fig 5: Phase Portrait: quiver and stream plots  
plt.rcParams["figure.figsize"] = (5,3)
xG = linspace(-1,1,10); yG = linspace(-5,5,5)
xx,yy = np.meshgrid(xG,yG)
xxDot = yy
yyDot = K1*xx+ K2*yy 

fig4, ax = plt.subplots(nrows=1, ncols=1)
txtx = 'x'; txty = 'v'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
fig4.tight_layout()
ax.plot(0,0,'ro',ms = 6)
ax.streamplot(xx,yy,xxDot,yyDot)
ax.plot(xP,yP,'b',lw = 2)

fig5, ax = plt.subplots(nrows=1, ncols=1)
txtx = 'x'; txty = 'v'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
fig5.tight_layout()
ax.quiver(xx,yy,xxDot,yyDot)
ax.plot(0,0,'ro',ms = 6)
ax.plot(xP,yP,'b',lw = 2)

#%%  Fig 6 Fig 7: LINEAR ALGEBRA SOLUTION
A = array([[0,1],[K1,K2]])
L, F = eig(A)
F0 = [[F[0,0],F[1,0]]]
F1 = [[F[0,1],F[1,1]]]
X = np.transpose(u0)
c = np.linalg.solve(F,X)

xA = c[0]*exp(L[0]*t)*F[0,0] + c[1]*exp(L[1]*t)*F[0,1]
vA = c[0]*exp(L[0]*t)*F[1,0] + c[1]*exp(L[1]*t)*F[1,1]

fig6, ax = plt.subplots(nrows=1, ncols=1)
txtx = 't'; txty = 'x'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.grid()
xP = t; yP = xA 
ax.plot(xP, yP,'b',lw = 2)
fig6.tight_layout()

fig7, ax = plt.subplots(nrows=1, ncols=1)
txtx = 't'; txty = 'v'
ax.set_xlabel(txtx);ax.set_ylabel(txty)
ax.grid()
xP = t; yP = v
ax.plot(xP,yP,'b',lw = 2)
fig7.tight_layout()

#%% CONSOLE OUTPUT
print(' ')
print('matrix A = ', np.round(A,2) )
print('eigenvalues = ', np.round(L,4))
print('eigenvector matrix = ', np.round(F,4))
print('eigenvector F0 = ', np.round(F0,4)) 
print('eigenvector F1 = ', np.round(F1,4) )

#%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')
fig6.savefig('a6.png')
fig7.savefig('a7.png')


'''
#