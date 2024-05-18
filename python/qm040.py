# -*- coding: utf-8 -*-
"""
qm040.py            May 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm040.pdf

QUANTUM MECHANICS
   Solution of the time independent Schrodinger equation
   for an electron bound in a potential well by finding the eigenvalues
   and eigenvectors
   

"""
# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt
import scipy.special 
import scipy.special as sps
from scipy import *
import scipy.linalg as la
import math
import cmath
import os
import matplotlib.pyplot as plt
from pylab import *
from math import *
from numpy import linspace,matrix,array
from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve
from scipy.integrate import odeint, solve_ivp, simps
from scipy.sparse.linalg import eigsh, eigs #Solves the Eigenvalue problem
from scipy.sparse import diags #Allows us to construct our matrices
import time

tStart = time.time()

#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
me = 9.10938356e-31            # Electron mass [kg]
se = e                         # Energy scaling factor   J <---> ev 
sx = 1e-9                      # Length scaling factor   m <---> nm


#%%  INPUTS
N = 519                #  grid size

xMin = -0.2*sx         #  default = -0.2 nm
xMax =  0.2*sx         #  default = +0.2 nm
U0 = -1000*se          #  Depth of well: default = -1000 eV
w = 0.1*sx             #  Width of well: default 0.1 nm
M = 30                 # number of eigenvalues returned
# >>> Enter 1,2,3,4,5,6 eigenstate number for expectation calculations
n = 1     


#%% COMPUTATIONS
x = linspace(xMin,xMax,N)           # x grid  [m]
dx = x[2] - x[1]

Cse = -hbar**2/(2*me) 
# Potential energy function  [J]
U = zeros(N)                
U[x>-w/2] = U0 
U[x>w/2] = 0                               
UM = diag(U)  
               
# AM (second derivative), KM (kinetic energy), HM (Hamiltonian) matrices
off = ones(N-1)
AM = (-2*np.eye(N) + np.diag(off,1) + np.diag(off,-1))/(dx**2)
KM = Cse*AM
HM = KM + UM

# Eigenvalues [J] and eigenfunctions (eigenvectors)
ev, ef = eigsh(HM, which="SM", k = M)


#%% OUTPUT: negative eigenvalues and normalize eigenfunctions
E = ev[ev<0]/se                 # negative eigenvalues [eV]
psi = zeros([N,len(E)]); psi2 = zeros([N,len(E)])
for c in range(len(E)):
    psi[:,c] = ef[:,c]
    psi2[:,c] = psi[:,c]**2
    area = simps(psi2[:,c],x)
    psi[:,c] = psi[:,c]/sqrt(area)

probD = psi**2    # probability density [1/m]
 
    
#%% GRAPHICS

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7)
fig1, axes = plt.subplots(nrows=3, ncols=2)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
                    right = 0.90, hspace = 0.36,wspace=0.40)

def graph(R,C,y,n,En,s):
    axes[R,C].xaxis.grid()
    axes[R,C].yaxis.grid()
    axes[R,C].set_title('n = %2.0f ' % n + '  E = %2.2f ev' % En, fontsize = 10)
    if s ==1:
       axes[R,C].set_ylabel('$\psi$  [a.u]',color = 'black',fontsize = 12)
       axes[R,C].set_ylim([-1.1, 1.1])
       axes[R,C].set_yticks([-1,-0.5,0,0.5,1])
    
    if s ==2:
       axes[R,C].set_ylabel('$|\psi|^2$ [1/m]',color = 'black',fontsize = 12)   
    
    x1 = -w/(2*sx); x2 = -x1
    axes[R,C].plot(x/sx,y, 'blue')
    axes[R,C].plot([x1,x1],[min(y),max(y)],'r')
    axes[R,C].plot([ x2, x2],[min(y),max(y)],'r')
    return    

# graph(R,c,psi,mode,E )

# Wavefunction plots
psiMax = amax(psi)
graph(0,0,psi[:,0]/psiMax, 1, E[0],1)
graph(0,1,psi[:,1]/psiMax, 2, E[1],1)
graph(1,0,psi[:,2]/psiMax, 3, E[2],1)
graph(1,1,psi[:,3]/psiMax, 4, E[3],1)
graph(2,0,psi[:,4]/psiMax, 5, E[4],1)
axes[2,0].set_xlabel('x  [nm]',color = 'black')
graph(2,1,psi[:,5]/psiMax, 6, E[5],1)
axes[2,1].set_xlabel('x  [nm]',color = 'black')

fig1.savefig('a1.png')


plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7)
fig1, axes = plt.subplots(nrows=3, ncols=2)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
                    right = 0.90, hspace = 0.36,wspace=0.40)

# graph(R,c,psi,mode,E )
graph(0,0,probD[:,0], 1, E[0],2)
graph(0,1,probD[:,1], 2, E[1],2)
graph(1,0,probD[:,2], 3, E[2],2)
graph(1,1,probD[:,3], 4, E[3],2)
graph(2,0,probD[:,4], 5, E[4],2)
axes[2,0].set_xlabel('x  [nm]',color = 'black')
graph(2,1,psi[:,5], 6, E[5],2)
axes[2,1].set_xlabel('x  [nm]',color = 'black')

fig1.savefig('a2.png')


#%% EXPECTATION VALUE CALCULATIONS

def firstDer(N,dx):
    v  = ones(N-1)
    M1 = diag(-v,-1)
    M2 = diag(v,1)
    M = M1+M2
    M[0,0] = -2; M[0,1] = 2; M[N-1,N-2] = -2; M[N-1,N-1] = 2
    MF = M/(2*dx) 
    return MF

def secondDer(N,dx):
     v = -2*ones(N)
     M1 = np.diag(v)
     v = np.ones(N-1)
     M2 = np.diag(v,1)
     M3 = np.diag(v,-1)
     M = M1+M2+M3
     M[0,0] = 1; M[0,1] = -2; M[0,2] = 1
     M[N-1,N-3] = 1; M[N-1,N-2] = -2; M[N-1,N-1]=1
     MS = M/(dx**2) 
     return MS

y = psi[:,n]      # eigenfunction n
# Probability
fn = y**2
prob = simps(fn,x)
# Position and its uncertainty
fn = y*x*y
x_avg = simps(fn,x)                 # [nm]
fn = y*x**2*y
x2_avg = simps(fn,x)                # [nm*nm]
dX = sqrt(x2_avg - x_avg**2)        #  [m]

# Momentum and its uncertainty
y2 = firstDer(N,dx)@y
fn = y*y2
p_avg = -1j*hbar*simps(fn,x)       # [N.s]
y2 = secondDer(N,dx)@y     # Second derivative matrix x [m]
fn = y*y2                  # Second derivative of function y
p2_avg = -hbar**2*simps(fn,x)   # [N^2.S^2]
dP = sqrt(p2_avg - imag(p_avg)**2)    # [N.s]

# Heisenberg Uncertainty Principle
HUP = 2*abs(dX*dP/hbar)

# Potential energy  [ev]
fn = y*U*y
U_avg = simps(fn,x)/se
# Kinetic energy [ev]
K_avg = p2_avg/(2*me)/se
# Total energy [ev]
E_avg = K_avg + U_avg


#%%  CONSOLE OUTPUT
print('  ')
print('grid point N = %2.0f' %N + '   eigenvalues returned M = %2.0f' %M)
s1 = xMin/sx; s2 = xMax/sx; print('xMin = %2.2f nm' % s1 + '   xMax = %2.2f nm' % s2)
s = w/sx; print('well width w = %2.2f nm' %s)  
s = U0/se; print('well depth U0 = %2.0f ev' %s)  
print(' ')
print('Energy eigenvalues [ev]')
for c in range(6):
    s = c+1; print('   E%0.0f' %s + ' = %2.3f' %E[c])
print(' ')  
print('Eigenstate n = %2.0f' % n)   
print('  Expectation values and Uncertainty Principle') 
print('    <x> = %2.2f m' %x_avg + '   deltax dX = %2.2e m' % dX)
s= real(p_avg); print('    <p> = %2.2f N.s' % s + '   deltax dP = %2.2e m' % dP)
print('    HUP = %2.2f ' % HUP + ' > 1') 
print('  ')
print('  Eigenstate energies' )
print('    En = %2.2f eV'  %E[n] )
print('    <E> = %2.2f' % E_avg + ' <K> = % 2.2f' % K_avg + ' <U> = % 2.2f' % U_avg)
s = K_avg + U_avg; print('    <K> + <U> = %2.2f' % s )


#%%  ENERGY PLOTS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)
fig, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('U  [ ev ]',color= 'black')
ax.set_xlabel('x  [nm]',color = 'black')
ax.plot(x/sx,U/se,'b',lw = 2)
xP = [xMin/sx,xMax/sx]
for c in range(6):
   yP = [E[c],E[c]]
   ax.plot(xP,yP,'r',lw = 2)

fig.tight_layout()

fig.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)

    

 
