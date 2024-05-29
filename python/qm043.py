# -*- coding: utf-8 -*-
"""
qm043.py    May 2024

QUANTUM MECHANICS
   Particle in a box
   Potential well with infinite depth

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation 
    
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm043.pdf


"""

#%%  IMPORTS
import numpy as np
from scipy.sparse import diags #Allows us to construct our matrices
from scipy.sparse.linalg import eigsh #Solves the Eigenvalue problem
import matplotlib.pyplot as plt
from scipy.linalg import *
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, real, imag
import time
from matplotlib.animation import FuncAnimation, PillowWriter 
from scipy.integrate import simps

tStart = time.time()



#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [ J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
me = 9.10938356e-31            # Electron mass [kg]
se = e                         # Energy scaling factor    J <---> eV
sx = 1e-9                      # Length scaling factor    m <---> m

#%% INPUTS: eigenstates n1 and n2 (quantum numbers n) / width of box
N = 10                # Number of eigenstates
nx = 599              # x grid points
n1 = 1            
n2 = 1
L = 0.1 * sx            # width of box  [m]

# ENERGY EIGENVALUES
n = np.arange(1,N+1)
E = n**2*(h**2/(8*me*L**2)) / se       #    [eV]

plt.rcParams["figure.figsize"] = (5,4)
plt.rcParams['font.size'] = 12
fig, ax = plt.subplots(1)
ax.set_ylabel('E  [ eV ]',color= 'black')
ax.set_xlabel('n',color = 'black')
ax.set_xticks(np.arange(1,11,1))
fig.tight_layout()
ax.plot(n,E,'bo')
#fig.savefig('a1.png')

for m in np.arange(1,N+1,1):
     # s = E[m-1]-1000; print('E%0.0f = ' % m + '%.3f eV' % s)
      s = E[m-1]; print('E%0.0f = ' % m + '%.3f eV' % s)
      
# NORMALIZED ENERGY EIGENVECTORS & PROABILITY DENSITY: normalized    
x = linspace(0,L, nx)
psi = zeros([N,nx]); probD = zeros([N,nx])
for m in range(N-1):
      psi[m+1,:] = sqrt(2/L) * sin((m+1)*pi*x/L)
      probD[m+1,:] = psi[m+1,:]**2

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
     #  axes[R,C].set_ylim([-1.1, 1.1])
     #  axes[R,C].set_yticks([-1,-0.5,0,0.5,1])
    
    if s ==2:
       axes[R,C].set_ylabel('$|\psi|^2$ [1/m]',color = 'black',fontsize = 12)   
    
    axes[R,C].plot(x/sx,y, 'blue')
    
    return    

# graph(R,c,psi,mode,E )
#%%
# Wavefunction plots
psiMax = np.amax(psi)
graph(0,0,psi[1,:]/psiMax, 1, E[0],1)
graph(0,1,psi[2,:]/psiMax, 2, E[1],1)
graph(1,0,psi[3,:]/psiMax, 3, E[2],1)
graph(1,1,psi[4,:]/psiMax, 4, E[3],1)
graph(2,0,psi[5,:]/psiMax, 5, E[4],1)
axes[2,0].set_xlabel('x  [nm]',color = 'black')
graph(2,1,psi[6,:]/psiMax, 6, E[5],1)
axes[2,1].set_xlabel('x  [nm]',color = 'black')

fig1.savefig('a1.png')


# Wavefunction plots
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7)
fig1, axes = plt.subplots(nrows=3, ncols=2)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
                    right = 0.90, hspace = 0.36,wspace=0.40)

graph(0,0,probD[1,:], 1, E[0],2)
graph(0,1,probD[2,:], 2, E[1],2)
graph(1,0,probD[3,:], 3, E[2],2)
graph(1,1,probD[4,:], 4, E[3],2)

graph(2,0,probD[5,:], 5, E[4],2)
axes[2,0].set_xlabel('x  [nm]',color = 'black')
graph(2,1,probD[6,:], 6, E[5],2)
axes[2,1].set_xlabel('x  [nm]',color = 'black')

fig1.savefig('a2.png')

#%%  ENERGY PLOTS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)
fig, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('E  [ ev ]',color= 'black')
ax.set_xlabel('x  [nm]',color = 'black')
ax.set_ylim([-1000,400])
xP = [0,L/sx]
for c in range(6):
   yP = [E[c]-1000,E[c]-1000]
   ax.plot(xP,yP,'r',lw = 2)

xP = [0,0]; yP = [-1000,400]
ax.plot(xP,yP,'b',lw = 4)
xP = [L/sx,L/sx]; yP = [-1000,400]
ax.plot(xP,yP,'b',lw = 4)
xP = [0,L/sx]; yP = [-1000,-1000]
ax.plot(xP,yP,'b',lw = 4)

fig.tight_layout()

fig.savefig('a3.png')