# -*- coding: utf-8 -*-
"""
qmCB.py    July 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm050.pdf

QUANTUM MECHANICS     
    Double well potential: COVALENT BONDING
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
from matplotlib.animation import FuncAnimation, PillowWriter 
import time

tStart = time.time()


#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
cL = 2.99792458e8              # speed of light
me = 9.10938356e-31            # electron mass [kg]
eps0 = 8.8541878188e-12        # episilon-zero  [F/m]
se = e                         # Energy scaling factor   J <---> ev 
sx = 1e-9                      # Length scaling factor   m <---> nm
ke = 1/(4*pi*eps0)             # Coulomb constant


#%%  INPUTS
xS = 1e-9
xMax = 5                   # X range limits  [nm] 
xMin = -5            
U0 = -100
N = 999                    # x grid points
num = 50                   # number of eigenvalues returned  


# SETUP: potential energy function
xs = linspace(xMin,xMax, N)         #  x grid   [nm]
#xe[abs(xe) < 0.005] = 0.005
x = xs*sx                            #  x grid   [m]

dx = x[2] - x[1]                        #  x grid spacing   [m]
U1 = -ke*e**2/abs(x+1e-18) 
U2 = -ke*e**2/abs(x-xS)                       #  potential energy  [J]
U = U1+U2
U = U2
#U[x>3e-9] = -2e-10
Us = U/se                               # potential energy  [eV] 

# FIG 1:  potential well plots
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,2)
fig1, ax = plt.subplots(1)


ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('U  [ eV ]',color= 'black')
ax.set_xlabel('x  [nm]',color = 'black')
ax.set_ylim([-50,10])
#ax.set_xlim([0,0.5])
ax.plot(xs,Us,'b',lw = 2)

fig1.tight_layout()


#%% EIGENVALUES AND EIGENFUNCTIONS
Cse = -hbar**2/(2*me)                   
UM = diag(U)  
               
# AM (second derivative), KM (kinetic energy), HM (Hamiltonian) matrices
off = ones(N-1)
AM = (-2*np.eye(N) + np.diag(off,1) + np.diag(off,-1))/(dx**2)
KM = Cse*AM
HM = KM + UM

# Eigenvalues [J] and eigenfunctions (eigenvectors)
ev, ef = eigsh(HM, which="SM", k = num)


#%% OUTPUT: negative eigenvalues and normalize eigenfunctions
E = ev[ev<0]/se                 # negative eigenvalues [eV]
psi = zeros([N,len(E)]); psi2 = zeros([N,len(E)])
for c in range(len(E)):
    psi[:,c] = ef[:,c]
    psi2[:,c] = psi[:,c]**2
    area = simps(psi2[:,c],x)
    psi[:,c] = psi[:,c]/sqrt(area)

probD = psi**2    # probability density [1/m]

lenE = len(E)

print(E)

#%% CONSOLE OUTPUT
print('  ')
# print('grid point N = %2.0f' % N + '   eigenvalues returned M = %2.0f' % num)
# # print('boundary width    xb = %2.2f nm' % xb )
# # print('well width        xw = %2.2f nm' % xw )
# # print('well separation   xs = %2.2f nm' % xs )
# # s = x0/sx; print('Equilibrium bond length x0 = %2.3f nm' %s)  
# # s = U0/se; print('well depth U0 = %2.3f ev' %s) 
# # print('spring constant k = %2.3f N/m' % k) 
# # print(' ')
print('Energies [eV]')
print('State n     E  [eV]')
for q in range(len(E)):
    s = q+1; print(' %0.0f' % s + '        %2.5f' % E[q] )
print(' ')  


#%% FIGS 2 & 3 eigenfunctions and probability distribution

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7.2)
fig2, axes = plt.subplots(nrows=3, ncols=2)
#fig2.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
#                     right = 0.90, hspace = 0.36,wspace=0.40)

def graph(R,C,y,n,En,s):
     axes[R,C].xaxis.grid()
     axes[R,C].yaxis.grid()
     q = n+1
     axes[R,C].set_title('n = %2.0f ' % q + '  E = %2.5f ev' % En, fontsize = 10)
     if s ==1:
         axes[R,C].set_ylabel('$\psi$  [a.u]',color = 'black',fontsize = 12)
         axes[R,C].set_ylim([-1.1, 1.1])
         axes[R,C].set_yticks([-1,-0.5,0,0.5,1])  
         axes[R,C].plot(xs,y, 'blue')   
      
    
     if s ==2:
         axes[R,C].set_ylabel('$|\psi|^2$ [1/m]',color = 'black',fontsize = 12)  
         axes[R,C].plot(xs,y, 'blue')
     return    

# graph(R,c,psi,mode,E )

# # FIG 2: Wavefunction 
psiMax = amax(psi)
graph(0,0,abs(psi[:,0]/psiMax), 0, E[0],1)
graph(0,1,psi[:,1]/psiMax, 1, E[1],1)
graph(1,0,psi[:,2]/psiMax, 2, E[2],1)
graph(1,1,psi[:,3]/psiMax, 3, E[3],1)
graph(2,0,psi[:,4]/psiMax, 4, E[4],1)
graph(2,1,psi[:,5]/psiMax, 5, E[5],1)

axes[2,0].set_xlabel('x  [nm]',color = 'black')
axes[2,1].set_xlabel('x  [nm]',color = 'black')
# #fig2.suptitle(' $x_S$ = %0.3f  nm' %xs )
fig2.tight_layout()



# # FIG 3:  probability density
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7)
fig3, axes = plt.subplots(nrows=3, ncols=2)
fig3.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
                    right = 0.90, hspace = 0.36,wspace=0.40)

# graph(R,c,psi,mode,E )
graph(0,0,probD[:,0], 0, E[0],2)
graph(0,1,probD[:,1], 1, E[1],2)
graph(1,0,probD[:,2], 2, E[2],2)
graph(1,1,probD[:,3], 3, E[3],2)
graph(2,0,probD[:,4], 4, E[4],2)
graph(2,1,probD[:,5], 5, E[5],2)

# #axes[2,0].set_xticks(arange(0,L+0.05,0.1))  
# #axes[2,1].set_xticks(arange(-0,L+0.05,0.1))  

axes[2,0].set_xlabel('x  [nm]',color = 'black')
axes[2,1].set_xlabel('x  [nm]',color = 'black')
#fig3.suptitle(' $x_S$ = %0.3f  nm' %xs )
fig3.tight_layout()


#%%  FIG 4: ENERGY PLOTS
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,2)
fig4, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('energy  [ eV ]',color= 'black')
ax.set_xlabel('x  [nm]',color = 'black')
ax.set_ylim([-20,2])

xP = [xMin,xMax]
for c in range(lenE):
    yP = [E[c],E[c]]
    col = [1,0,0]
    q = c%2
    if q == 1:
        col = [0,0,0]
    ax.plot(xP,yP,lw = 2,color = col)

ax.plot(xs,Us,'b',lw = 3)
fig4.tight_layout()


#%% ENERGY as a function of separation distance
XS = array([0,0.7,0.1])
# Symmetric states
ES = array([-48.77711,-48.78364,-48.79192])    
# Antisymmetric states      
EA = array([-48.77711,-48.78364,-48.77806])           


#%%  SAVE FIGURES
# fig1.savefig('a1.png')
# fig2.savefig('a2.png')
# fig3.savefig('a3.png')
# fig4.savefig('a4.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)

