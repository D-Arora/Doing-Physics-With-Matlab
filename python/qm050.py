# -*- coding: utf-8 -*-
"""
qm050.py    June 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm050.pdf

QUANTUM MECHANICS     HCl molecule
   Solution of the time independent Schrodinger equation
   for an electron bound in a potential well by finding the eigenvalues
   and eigenvectors
   Parabolic potential well and Morse potential
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
mu = 1.62747e-27               # reduced mass HCl molecule [kg]
se = e                         # Energy scaling factor   J <---> ev 
sx = 1e-9                      # Length scaling factor   m <---> nm

#%%  INPUTS
N = 801                     # x grid points
xMin = 0; xMax = 0.3*sx     # x range [m]
L = 0.3                     # Limit for X axis in plots
x0 = 0.127*sx               # Equilibrium bond length HCl  [m]
U0 = -4.57*se               # Well depth [J]
k = 480;                    # HCL molecule: spring constant [N/m]

num = 50                    # number of eigenvalues returned  

m = 4; n = 3                # quantum numbers 0,1,2,3, ...(M > n)

am = 1; an = 0.5            # Eigenstate coeffs for compound state
s = sqrt(am**2 + an**2)
am = am/s; an = an/s        # normalized values

#%% SETUP: potential energy function
x = linspace(xMin, xMax, N)          # x grid
dx = x[2] - x[1] 
Uh = 0.5*k*(x-x0)**2 + U0            # Parabolic potential well  [J]
#
S = sqrt(2*abs(U0)/k)                # Width of Morse potential [m]
Um = U0*(1 - (exp((x0 - x)/S)-1)**2)

# FIG 1:  potential well plots
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('U  [ eV ]',color= 'black')
ax.set_xlabel('bond length x  [nm]',color = 'black')
ax.set_ylim([-6.1,6.1])
ax.set_xlim([0,0.5])
ax.plot(x/sx,Uh/se,'b',lw = 2)
ax.plot(x/sx,Um/se,'r',lw = 2)
fig1.tight_layout()

fig1.savefig('a1.png')


#%% EIGENVALUES AND EIGENFUNCTIONS
Cse = -hbar**2/(2*mu) 
U = Uh                           
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


# Theoretical eignvalues ET parabolic well & calculated binding energies EB
lenE = len(E)
ET = zeros(lenE)
omega = sqrt(k/mu)
for q in range(lenE):
    ET[q] = (q+0.5)*hbar*omega/se       #    [eV]
Ew = E - U0/se
EB = -E
# Spacing between adjacent energy levels
dE = zeros(lenE)
for q in range(len(E)-1):
    dE[q+1] = E[q+1]-E[q]
dE[0]= 0


#%% CONSOLE OUTPUT
print('  ')
print('grid point N = %2.0f' % N + '   eigenvalues returned M = %2.0f' % num)
s1 = xMin/sx; s2 = xMax/sx; print('xMin = %2.2f nm' % s1 + '   xMax = %2.2f nm' % s2)
s = x0/sx; print('Equilibrium bond length x0 = %2.3f nm' %s)  
s = U0/se; print('well depth U0 = %2.3f ev' %s) 
print('spring constant k = %2.3f N/m' % k) 
print(' ')
print('Energies [eV]')
print('State n    Ew       ET      dE       E       EB')
for q in range(len(E)):
    print(' %0.0f' % q + '        %2.3f' % Ew[q] + '    %2.3f' %ET[q] 
          + '   %2.3f' % dE[q] + '   %2.3f' %E[q] +  '   %2.3f'  %EB[q])
print(' ')  


#%% FIGS 2 & 3 eigenfunctions and probability distribution

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7.2)
fig2, axes = plt.subplots(nrows=3, ncols=2)
fig2.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
                    right = 0.90, hspace = 0.36,wspace=0.40)

def graph(R,C,y,n,En,s):
    axes[R,C].xaxis.grid()
    axes[R,C].yaxis.grid()
    axes[R,C].set_xticks(arange(0,L+0.05,0.1)) 
    axes[R,C].set_title('n = %2.0f ' % n + '  E = %2.2f ev' % En, fontsize = 10)
    if s ==1:
       axes[R,C].set_ylabel('$\psi$  [a.u]',color = 'black',fontsize = 12)
       axes[R,C].set_ylim([-1.1, 1.1])
       axes[R,C].set_yticks([-1,-0.5,0,0.5,1])  
       
       axes[R,C].plot(x/sx,y, 'blue')
       
    #   axes[R,C].plot(x/sx,-0.5+U/(se*10),'r')   
      
    if s ==2:
       axes[R,C].set_ylabel('$|\psi|^2$ [1/m]',color = 'black',fontsize = 12)  
       
       axes[R,C].plot(x/sx,y, 'blue')
       
      # axes[R,C].plot(x/sx,0.5+1e9*U/(se),'r')   
       
    return    

# graph(R,c,psi,mode,E )

# FIG 2: Wavefunction 
psiMax = amax(psi)
graph(0,0,abs(psi[:,0]/psiMax), 0, E[0],1)
graph(0,1,psi[:,1]/psiMax, 1, E[1],1)
graph(1,0,psi[:,2]/psiMax, 2, E[2],1)
graph(1,1,psi[:,3]/psiMax, 3, E[3],1)
graph(2,0,psi[:,4]/psiMax, 4, E[4],1)
graph(2,1,psi[:,5]/psiMax, 5, E[5],1)

axes[2,0].set_xlabel('x  [nm]',color = 'black')
axes[2,1].set_xlabel('x  [nm]',color = 'black')
fig1.tight_layout()

fig2.savefig('a2.png')

# FIG 3:  probability density
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

axes[2,0].set_xticks(arange(0,L+0.05,0.1))  
axes[2,1].set_xticks(arange(-0,L+0.05,0.1))  

axes[2,0].set_xlabel('x  [nm]',color = 'black')
axes[2,1].set_xlabel('x  [nm]',color = 'black')
fig3.tight_layout()
fig3.savefig('a3.png')


#%%  FIG 4: ENERGY PLOTS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('energy  [ eV ]',color= 'black')
ax.set_xlabel('x  [nm]',color = 'black')
ax.set_ylim([-5,0])

xP = [xMin/sx,xMax/sx]
for c in range(lenE):
   yP = [E[c],E[c]]
   ax.plot(xP,yP,'r',lw = 2)

ax.plot(x/sx,U/se,'b',lw = 3)

fig4.tight_layout()

fig4.savefig('a4.png')


#%% Eigensates m and n
# FIG 5: Eigenfunction plot for state n

def fig56(yP,q):
   # plt.rcParams["figure.figsize"] = (6,3)
   # fig, ax = plt.subplots(1)
   # fig.subplots_adjust(top = 0.90, bottom = 0.20, left = 0.18,\
   #                  right = 0.80, hspace = 0.20,wspace=0.2)

   ax.set_xticks(arange(-0.5,0.6,0.25)) 
   ax.xaxis.grid()
   ax.yaxis.grid()
   ax.set_ylabel('$\psi$  [a.u.]',color= 'blue')
   ax.set_xlabel('x  [nm]',color = 'black')
   ax.tick_params(axis='y', labelcolor="blue")
   SS = 1
   if yP[50] < 0:
      SS = -1
   xP = x/sx; yP = SS*yP
   ax.plot(xP,yP,'b',lw = 2)
   s = q; ax.set_title('E%0.0f'%q + '= %2.3f eV' %E[q])
   ax.set_yticks(arange(-1,1.2,0.5)) 
   ax.set_xticks(arange(0,L+0.05,0.05))  

   ax2 = ax.twinx()
   ax2.set_ylabel('U & K [eV]',color= 'k')
   ax2.tick_params(axis='y', labelcolor = 'k')
   ax2.plot(x/sx,U/se,'r',lw = 1)
   ax2.plot(x/sx,E[n]-U/se,'g',lw = 1)
   ax2.set_ylim([-20,20])
   ax2.set_yticks(np.arange(-20,21,5)) 

q = m; yP = psi[:,m]/amax(psi)
plt.rcParams["figure.figsize"] = (6,3)
fig5, ax = plt.subplots(1)
fig5.subplots_adjust(top = 0.90, bottom = 0.20, left = 0.18,\
                 right = 0.80, hspace = 0.20,wspace=0.2)

fig56(yP,q)

fig5.savefig('a5.png')



q = n; yP = psi[:,n]/amax(psi)
plt.rcParams["figure.figsize"] = (6,3)
fig6, ax = plt.subplots(1)
fig6.subplots_adjust(top = 0.90, bottom = 0.20, left = 0.18,\
                 right = 0.80, hspace = 0.20,wspace=0.2)
fig56(yP,q)

fig6.savefig('a6.png')

def fig78(q):
    
    ax.set_xticks(arange(-0.5,0.6,0.25)) 
    ax.xaxis.grid()
    ax.yaxis.grid()
    ax.set_ylabel('$|\psi|^2$  [1/m]',color= 'blue')
    ax.set_xlabel('x  [nm]',color = 'black')
    ax.tick_params(axis='y', labelcolor="blue")
    ax.plot(x/sx,probD[:,q],'b',lw = 2)
    s = q; ax.set_title('E%0.0f'%q + '= %2.3f eV' %E[q])
    ax.set_xticks(arange(0,L+0.05,0.05))  
    
    ax2 = ax.twinx()
    ax2.set_ylabel('U & K [eV]',color= 'k')
    ax2.tick_params(axis='y', labelcolor = 'k')
    ax2.plot(x/sx,U/se,'r',lw = 1)
    ax2.plot(x/sx,E[q]-U/se,'g',lw = 1)
    ax2.set_ylim([-20,20])
    ax2.set_yticks(np.arange(-20,21,5))
        
q = m
plt.rcParams["figure.figsize"] = (6,3)

fig7, ax = plt.subplots(1)
fig7.subplots_adjust(top = 0.90, bottom = 0.20, left = 0.18,\
                    right = 0.80, hspace = 0.20,wspace=0.2)
fig78(q)

fig7.savefig('a7.png')

q = n
plt.rcParams["figure.figsize"] = (6,3)

fig8, ax = plt.subplots(1)
fig8.subplots_adjust(top = 0.90, bottom = 0.20, left = 0.18,\
                    right = 0.80, hspace = 0.20,wspace=0.2)
fig78(q)

fig8.savefig('a8.png')    


#%%  COMPOUND STATE  m and n
PSI = am*psi[:,m] + an*psi[:,n]
PROBD = PSI**2
PROB = simps(PROBD,x)        # Check probability = 1

# Expection value <x>
fn = PSI*x*PSI
xavg = simps(fn,x)/sx

wm = Ew[m]*se/hbar; wn = E[n]*se/hbar
wmn = abs(wm-wn)
Tmn = 2*pi/wmn
fmn = 1/Tmn
Lmn = cL/fmn

y = am*psi[:,m]*exp(-1j*wm*Tmn/2) + an*psi[:,n]*exp(-1j*wn*Tmn/2)
yC = conj(y)
fn = real(yC*x*y)
xavg1 = simps(fn,x)/sx
xD1 = x0/sx - xavg
xD2 = x0/sx - xavg1
xD = abs(xD1)+abs(xD2)

# FIG 9:  Compound state
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig9, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('$\Psi$  [a.u.]',color= 'black',fontsize = 12)
ax.set_title('m = %0.0f' %m + ' am = %0.2f' % am + 
         '   n = %0.0f' % n + ' an = %0.2f' % an, fontsize = 12)
ax.set_xlabel('x  [nm]',color = 'black')
ax.set_xticks(arange(0,L+0.05,0.05))  
xP = x/sx; yP = PSI/max(PSI)
ax.plot(xP, yP, 'blue')
fig9.tight_layout()

fig9.savefig('a9.png')

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig10, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('$|\Psi|^2$  [1/m]',color= 'black',fontsize = 12)

ax.set_title('m = %0.0f' %m + ' am = %0.2f' % am + 
         '   n = %0.0f' % n + ' an = %0.2f' % an, fontsize = 12)
ax.set_xlim([0.06,0.18])
#ax.set_xticks(arange(0,L+0.05,0.05))  
xP = x/sx; yP = PROBD
ax.plot(xP, yP, 'blue',lw = 2)
yP = real(yC*y)
ax.plot(xP, yP, 'red',lw = 2)
xP = [x0/sx,x0/sx];
yP = [0,1.2*max(PROBD)]
ax.plot(xP, yP, 'm',lw=5) 
xP = [xavg,xavg] 
ax.plot(xP, yP, 'k',lw=2) 
xP = [xavg1,xavg1]
ax.plot(xP, yP, 'k',lw=2) 
                         
ax.set_xlabel('x  [nm]',color = 'black',fontsize = 12)
fig10.tight_layout()

fig10.savefig('a10.png')

#%% CONSOLE: Compound state
print('  ')
print('COMPOUND STATE')
print('m = %0.0f' % m + ' am = %2.3f' % am + '   n = %0.0f' %n 
      + ' an = %2.3f' % an)
s = E[m] - E[n]; print('dE = %2.3f  eV' %s)      
print('osc. frequency omega = %2.3e  rad/s' %wmn) 
print('osc. frequency f = %2.3e  Hz' %fmn) 
print('osc. period T = %2.3e  s' %Tmn)
s = Lmn/sx;print('photon wavelength lambda = %2.3f  nm' %s)
s = x0/sx; print('equilibrium bond length = %2.3f  nm' %s)
print('el. dipole separation = %2.7f  nm' %xD)

#%% SCATTER PLOT FOR PROBILITY
q = n
y = probD[:,q]/amax(probD)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,1)
fig11, ax = plt.subplots(1)
fig11.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.18,\
                    right = 0.80, hspace = 0.20,wspace=0.2)

x1 = xMin/sx; x2 = xMax/sx
ax.set_xlim([x1,x2])
ax.set_ylim([0,1])
ax.set_axis_off()

random.seed()
num = 20000
for c in range(num):
    xP = x1+x2*random.random()        # [mm]
    yP = random.random()
    qR = random.random()
    qP = y[x/sx>xP]
    q0 = qP[0] 
    if qR < q0:
       ax.plot(xP,yP,'bo', ms = 1)
fig11.savefig('a11.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)

    