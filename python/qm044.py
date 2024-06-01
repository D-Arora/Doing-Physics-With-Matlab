# -*- coding: utf-8 -*-
"""
qm044.py            June 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm044.pdf

QUANTUM MECHANICS
   Solution of the time independent Schrodinger equation
   for an electron bound in a potential well by finding the eigenvalues
   and eigenvectors
   Sloping floor finite potential well

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
U1 = -1200*se          #  Depth of well: default = -1200 eV
U2 = -200*se
w = 0.2*sx             #  Width of well: default 0.2 nm
M = 30                 # number of eigenvalues returned
# >>> Enter 1,2,3,4,5,6 eigenstate number for expectation calculations and plots
n = 6     


#%% COMPUTATIONS
x = linspace(xMin,xMax,N)           # x grid  [m]
dx = x[2] - x[1]

Cse = -hbar**2/(2*me) 

# Potential energy function  [J]
U = zeros(N) 
x1 = -(w/2); x2 = (w/2)
m = (U2-U1)/(x2-x1); b = U1 - m*x1             
U[x>-w/2] = m*x[x>-w/2]  + b 
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
plt.rcParams["figure.figsize"] = (7,7.2)
fig1, axes = plt.subplots(nrows=3, ncols=2)
fig1.subplots_adjust(top = 0.95, bottom = 0.07, left = 0.120,\
                    right = 0.98, hspace = 0.60,wspace=0.40)

def graph(R,C,y,n,En,s):
    axes[R,C].xaxis.grid()
    axes[R,C].yaxis.grid()
    axes[R,C].set_title('n = %2.0f ' % n + '  E = %2.2f ev' % En, fontsize = 10)
    if s ==1:
       axes[R,C].set_ylabel('$\psi$  [a.u]',color = 'black',fontsize = 12)
       axes[R,C].set_ylim([-1.1, 1.1])
       axes[R,C].set_yticks([-1,-0.5,0,0.5,1])
       axes[R,C].set_xlabel('x  [nm]',color = 'black')
       #x1 = -w/(2*sx); x2 = -x1
       axes[R,C].plot(x/sx,U/abs(U1),'r',lw = 2)
       
       S = 1
       if y[50] < 0:
          S = -1
       xP = x/sx; yP = S*y
       axes[R,C].plot(xP,yP,'b',lw = 2)
       
    if s ==2:
       axes[R,C].set_ylabel('$|\psi|^2$ [1/m]',color = 'black',fontsize = 12)   
       axes[R,C].plot(x/sx,y, 'blue')
       axes[R,C].set_xlabel('x  [nm]',color = 'black')
       axes[R,C].set_xticks(arange(-0.2,0.3,0.1)) 
       
       axes[R,C].plot(x/sx,1.2e10 + 1.2e10*U/abs(U1),'r',lw = 2)
    return    

# graph(R,c,psi,mode,E )

# FIG 1: Wavefunction plots
s = zeros(2)
#s[0] = amax(psi); s[1] = abs(amin(psi)); s = max(s)
s[0] = max(psi[:,0]); s[1] = abs(min(psi[:,0])); s = max(s)
 
graph(0,0,psi[:,0]/s, 1, E[0],1)
graph(0,1,psi[:,1]/s, 2, E[1],1)
graph(1,0,psi[:,2]/s, 3, E[2],1)
graph(1,1,psi[:,3]/s, 4, E[3],1)
graph(2,0,psi[:,4]/s, 5, E[4],1)
graph(2,1,psi[:,5]/s, 6, E[5],1)

fig1.savefig('a1.png')

# FIG 2: Probability densities
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7)
fig1, axes = plt.subplots(nrows=3, ncols=2)
fig1.subplots_adjust(top = 0.95, bottom = 0.07, left = 0.120,\
                    right = 0.98, hspace = 0.54,wspace=0.40)
# graph(R,c,psi,mode,E )
graph(0,0,probD[:,0], 1, E[0],2)
graph(0,1,probD[:,1], 2, E[1],2)
graph(1,0,probD[:,2], 3, E[2],2)
graph(1,1,probD[:,3], 4, E[3],2)
graph(2,0,probD[:,4], 5, E[4],2)
graph(2,1,probD[:,5], 6, E[5],2)

# axes[2,0].set_xticks(arange(-0.2,0.3,0.1))  
# axes[2,1].set_xticks(arange(-0.2,0.3,0.1))  

# axes[2,0].set_xlabel('x  [nm]',color = 'black')
# axes[2,1].set_xlabel('x  [nm]',color = 'black')

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

y = psi[:,n-1]      # eigenfunction n
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
#s = U0/se; print('well depth U0 = %2.0f ev' %s)  
print(' ')
print('Energy eigenvalues [ev]')
for c in range(len(E)):
    s = c+1; print('   E%0.0f' %s + ' = %2.3f' %E[c])
print(' ')  
print('Eigenstate n = %2.0f' % n)   
print('  Expectation values and Uncertainty Principle') 
s = x_avg/sx; print('    <x> = %2.3f nm' %s) 
s = dX/sx;    print('    deltax dX = %2.3f nm' % s)
s = real(p_avg); print('    <p> = %2.2f N.s' % s + '   deltax dP = %2.2e m' % dP)
print('    HUP = %2.2f ' % HUP + ' > 1') 
print('  ')
print('  Eigenstate energies' )
print('    En = %2.2f eV'  %E[n-1] )
print('    <E> = %2.2f' % E_avg + ' <K> = % 2.2f' % K_avg + ' <U> = % 2.2f' % U_avg)
s = K_avg + U_avg; print('    <K> + <U> = %2.2f' % s )


#%%  FIG 3: ENERGY PLOTS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4.5)
fig, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('energy  [ ev ]',color= 'black')
ax.set_xlabel('x  [nm]',color = 'black')
ax.set_ylim([-1200,50])
ax.plot(x/sx,U/se,'b',lw = 2)

xP = [xMin/sx,xMax/sx]
for c in range(len(E)):
   yP = [E[c],E[c]]
   ax.plot(xP,yP,'r',lw = 2)

ax.plot(x/sx,U/se,'b',lw = 3)

fig.tight_layout()

fig.savefig('a3.png')


#%%   EIGENVALUE EXPLORATIONS   state n
y = probD[:,n-1]/amax(probD)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,1)
fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.18,\
                    right = 0.80, hspace = 0.20,wspace=0.2)

x1 = xMin/sx; x2 = xMax/sx
ax.set_xlim([x1,x2])
ax.set_ylim([0,1])
ax.set_axis_off()

random.seed()
num = 20000
for c in range(num):
    xP = x1+2*x2*random.random()        # [mm]
    yP = random.random()
    qR = random.random()
    qP = y[x/sx>xP]
    q0 = qP[0] 
    if qR < q0:
       ax.plot(xP,yP,'bo', ms = 1)
fig.savefig('a4.png')


#%% Eigenfunction plot for state n
plt.rcParams["figure.figsize"] = (6,3)

fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.18,\
                    right = 0.80, hspace = 0.20,wspace=0.2)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('$\psi$  [a.u.]',color= 'blue')
ax.set_xlabel('x  [nm]',color = 'black')
ax.tick_params(axis='y', labelcolor="blue")
ax.plot(x/sx,y,'b',lw = 2)

ax2 = ax.twinx()
ax2.set_ylabel('U & K [eV]',color= 'k')
ax2.set_yticks(np.arange(-1200,410,400))
ax2.tick_params(axis='y', labelcolor = 'k')
ax2.plot(x/sx,U/se,'r',lw = 1)
ax2.plot(x/sx,E[n-1]-U/se,'g',lw = 1)
fig.savefig('a5.png')


#%% Probability density plot for state n
plt.rcParams["figure.figsize"] = (6,3)

fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.18,\
                    right = 0.80, hspace = 0.20,wspace=0.2)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('$|\psi|^2$  [1/m]',color= 'blue')
ax.set_xlabel('x  [nm]',color = 'black')
ax.tick_params(axis='y', labelcolor="blue")
ax.plot(x/sx,probD[:,n-1],'b',lw = 2)

ax2 = ax.twinx()
ax2.set_ylabel('U & K [eV]',color= 'k')
ax2.set_yticks(np.arange(-1200,410,400))
ax2.tick_params(axis='y', labelcolor = 'k')
ax2.plot(x/sx,U/se,'r',lw = 1)
ax2.plot(x/sx,E[n-1]-U/se,'g',lw = 1)

fig.savefig('a6.png')

    
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)

    

 
