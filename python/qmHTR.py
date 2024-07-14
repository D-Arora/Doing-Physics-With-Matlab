# -*- coding: utf-8 -*-

"""
qmHTR.py    July 2024

Ian Cooper 
        matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm061.pdf

QUANTUM MECHANICS     HYDEOGEN ATOM: SELECTIN RULES   TRANSITION RATES
  
"""

# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, radians
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
from scipy.special import sph_harm

tStart = time.time()

#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
cL = 2.99792458e8              # speed of light
eps0 = 8.8541878188e-12        # episilon-zero  [F/m]
me = 9.109e-31
mp = 1.673e-27
mu = me*mp/(me+mp)             # reduced mass HCl molecule [kg]
se = e                         # Energy scaling factor   J <---> ev 
sx = 1e-9                      # Length scaling factor   m <---> nm
ke = 1/(4*pi*eps0)             # Coulomb constant


#%%  INPUTS
# INITIAL eigenstate n1, L1, mL1  
n1 = 7; L1 = 1; mL1 = 0
# FINAL eigenstate n2, L2, mL2   
n2 = 2; L2 = 0; mL2 = 0

rMax = 7*sx       # r range max [m]

Z = 1               # nuclear charge
#N = 255             # r grid points
num = 100           # number of eigenvalues returned  
rMin = 1e-18        # r range min  [m]

N = 999


print('INPUTS  ')
s1 = L2 - L1; s2 = mL2 - mL1
print('  L1 = %0.0f' %L1 + '  L2 = %0.0f' %L2 +
            ' DeltaL = %0.0f' %s1  + '   mL1 = %0.0f' %mL1 +      
            '   mL2 = %0.0f' %mL2 + r'  DeltamL = %0.0f' %s2)
s1 = rMax/sx; print('  rMax = %2.3f nm' %s1 + '   grid points N = %0.0f' %N)
print(' ')


#%%   SELECTION RULES: MAGNETIC QUANTUM NUMBER mL
# Arimuthal angle  [rad]
phi = linspace(0,2*pi,N)        
# Azimuthal wavefunction
PHI1 = 1/sqrt(2*pi) * exp(-1j*mL1*phi)
PHI2C = 1/sqrt(2*pi) * exp(1j*mL2*phi)
# Integrands
fAx = PHI2C*PHI1*cos(phi)
fAy = PHI2C*PHI1*sin(phi)
fAz = PHI2C*PHI1
IAx = simps(fAx,phi)
IAy = simps(fAy,phi)
IAz = simps(fAz,phi)

print('Azimuthal equations ')
s1 = IAx.real; s2 = IAx.imag 
print('IAx = %0.3f' %s1 + ' +  %0.3fj' %s2) 
s1 = IAy.real; s2 = IAy.imag 
print('IAy = %0.3f' %s1 + ' +  %0.3fj' %s2)            
s1 = IAz.real; s2 = IAz.imag 
print('IAz = %0.3f' %s1 + ' +  %0.3fj' %s2) 


#%%  FIG 1: Azimuthal functions
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=3)

def graph(C,fn,txt,txtT):
    ax[C].plot(phi/pi, fn,color = col,lw = 2)
    ax[C].fill_between(phi/pi,fn,color = col)
    ax[C].set_xlabel(r'$\phi / \pi$', fontsize = 12)
    ax[C].set_ylabel(txt, fontsize = 12)
    ax[C].grid(True)
    ax[C].set_ylim([-0.2,0.2]) 
    ax[C].set_xlim([0,2])
    ax[C].set_xticks(arange(0,2.2,0.5))
    ax[C].set_title(txtT,fontsize = 10)
      
txt = 'f$_{Ax}$'
s1 = IAx.real; s2 = IAx.imag; txtT = 'IAx = %0.2f' %s1 + ' + %0.2fj' %s2
fn = fAx.real; col = [0,0,1]
if abs(IAx.imag) > abs(IAx.real):
    fn = fAx.imag; col = [1,0,0]
if abs(IAx) < 0.1:
    col = [1,0,1]
graph(0,fn,txt,txtT)

txt = 'f$_{Ay}$'
s1 = IAy.real; s2 = IAy.imag; txtT = 'IAy = %0.2f' %s1 + ' + %0.2fj' %s2
fn = fAy.real; col = [0,0,1]
if abs(IAy.imag) > abs(IAy.real):
    fn = fAy.imag; col = [1,0,0]
if abs(IAy) < 0.1:
    col = [1,0,1]
graph(1,fn,txt,txtT)

txt = 'f$_{Az}$'
s1 = IAz.real; s2 = IAz.imag; txtT = 'IAz = %0.2f' %s1 + ' + %0.2fj' %s2
fn = fAz.real; col = [0,0,1]
if abs(IAz.imag) > abs(IAz.real):
    fn = fAz.imag; col = [1,0,0]
if abs(IAz) < 0.1:
    col = [1,0,1]    
graph(2,fn.real,txt,txtT)

s = mL2-mL1
plt.suptitle('  mL$_1$ = %0.0f' %mL1 +      
         '    mL$_2$ = %0.0f' %mL2 + r'    $\Delta$mL = %0.0f' %s)
fig1.tight_layout() 


#%%   SELECTION RULES: ORBITAL ANGUALR MOMENTUM QUANTUM NUMBER L
theta = linspace(0,pi,N)
# Associated Legendre functions
THETA1 = sph_harm(mL1, L1, 0, theta).real
THETA2 = sph_harm(mL2, L2, 0, theta).real
# Normalize wavefunctions
fn = THETA1**2*sin(theta)
area = simps(fn,theta)
THETA1 = THETA1/sqrt(area)

fn = THETA2**2*sin(theta)
area = simps(fn,theta)
THETA2 = THETA2/sqrt(area)

# Integrands
fSx = THETA1*THETA2*sin(theta)**2
fSy = THETA1*THETA2*sin(theta)**2
fSz = THETA1*THETA2*cos(theta)*sin(theta)
ISx = simps(fSx,theta)
ISy = simps(fSy,theta)
ISz = simps(fSz,theta)

print('  ')
print('Polar equations ')
s1 = ISx.real; s2 = ISx.imag 
print('ISx = %0.3f' %s1 + ' +  %0.3fj' %s2) 
s1 = ISy.real; s2 = ISy.imag 
print('ISy = %0.3f' %s1 + ' +  %0.3fj' %s2)            
s1 = ISz.real; s2 = ISz.imag 
print('ISz = %0.3f' %s1 + ' +  %0.3fj' %s2) 


#%%  FIG 2: Polar functions
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=3)

def graph2(C,fn,txt,txtT):
    ax[C].plot(theta/pi, fn,color = col,lw = 2)
    ax[C].fill_between(theta/pi,fn,color = col)
    ax[C].set_xlabel(r'$\theta / \pi$', fontsize = 12)
    ax[C].set_ylabel(txt, fontsize = 12)
    ax[C].grid(True)
   # ax[C].set_ylim([-0.2,0.2]) 
   # ax[C].set_yticks(arange(-0.6,0.6,0.5))
    ax[C].set_xlim([0,1])
    ax[C].set_xticks(arange(0,1.1,0.5))
    ax[C].set_title(txtT,fontsize = 10)
    
txt = 'f$_{Sx}$'
s1 = ISx.real; s2 = ISx.imag; txtT = 'ISx = %0.2f' %s1 + ' + %0.2fj' %s2
fn = fSx.real; col = [0,0,1]
if abs(ISx.imag) > abs(ISx.real):
    fn = fSx.imag; col = [1,0,0]
if abs(ISx) < 0.1:
    col = [1,0,1]
graph2(0,fn,txt,txtT)


txt = 'f$_{Sy}$'
s1 = ISy.real; s2 = ISy.imag; txtT = 'ISy = %0.2f' %s1 + ' + %0.2fj' %s2
fn = fSy.real; col = [0,0,1]
if abs(ISy.imag) > abs(ISy.real):
    fn = fSy.imag; col = [1,0,0]
if abs(ISy) < 0.1:
    col = [1,0,1]
graph2(1,fn,txt,txtT)

txt = 'f$_{Az}$'
s1 = ISz.real; s2 = ISz.imag; txtT = 'ISz = %0.2f' %s1 + ' + %0.2fj' %s2
fn = fSz.real; col = [0,0,1]
if abs(ISz.imag) > abs(ISz.real):
    fn = fSz.imag; col = [1,0,0]
if abs(ISz) < 0.1:
    col = [1,0,1]    
graph2(2,fn.real,txt,txtT)

s = mL2-mL1; sL = L2-L1
plt.suptitle('L$_1$ = %0.0f' %L1 + '  L$_2$ = %0.0f' %L2 +
             r'  $\Delta$L = %0.0f' %sL +
            '     mL$_1$ = %0.0f' %mL1 +      
            '    mL$_2$ = %0.0f' %mL2 + r'   $\Delta$mL = %0.0f' %s)
fig2.tight_layout() 


#%%   COMPONENTS FOR THE OVERLAP ANGULAR INTEGRALS 
IASx = IAx*ISx
IASy = IAy*ISy
IASz = IAz*ISz
# Console output
print('  ')
print('Angular equations')
s1 = IASx.real; s2 = IASx.imag 
print('IASx = %0.3f' %s1 + ' +  %0.3fj' %s2) 
s1 = IASy.real; s2 = IASy.imag 
print('IASy = %0.3f' %s1 + ' +  %0.3fj' %s2)            
s1 = IASz.real; s2 = IASz.imag 
print('IASz = %0.3f' %s1 + ' +  %0.3fj' %s2) 


#%% RADIAL EQUATIONS
if n1 <= L1:
    print('  ')
    print('EXIT: n>L not satisfied')
    sys.exit(0)

# Potential energy functions
r = linspace(rMin,rMax,N)     # radial distance [m]
dx = r[2] - r[1]
# Coulomb potential energy function   [J]
UC = -ke*Z*e**2/r                   
# Orbital kinetic energy function     [J]
UL1 = L1*(L1+1)*hbar**2/(2*mu*r**2)
UL2 = L2*(L2+1)*hbar**2/(2*mu*r**2)
# Effective potential energy function  [J]
U1 = UC + UL1; U2 = UC + UL2

#%% EIGENVALUES AND EIGENFUNCTIONS
Cse = -hbar**2/(2*mu)               # Schrodinger constant
UM1 = diag(U1); UM2 = diag(U2)                       # potential energy matrix
               
# AM (second derivative), KM (kinetic energy), HM (Hamiltonian) matrices
off = ones(N-1)
AM = (-2*np.eye(N) + np.diag(off,1) + np.diag(off,-1))/(dx**2)
KM = Cse*AM

# STATE 1
HM1 = KM + UM1
# Eigenvalues [J] and eigenfunctions (eigenvectors)
ev, ef = eigsh(HM1, which="SM", k = num)
E = real(ev[ev<0]/se)                 # negative eigenvalues [eV]
lenE = len(E)
EB = -E
N1 = n1-L1-1
EB1 = EB[N1]

psi = zeros([N,len(E)]); psi2 = zeros([N,len(E)])
for c in range(len(E)):
    psi[:,c] = ef[:,c]
    psi2[:,c] = psi[:,c]**2
    area = simps(psi2[:,c],r)
    psi[:,c] = psi[:,c]/sqrt(area)
probD1 = psi**2    # probability density [1/m]
g1 = psi[:,N1]

print('  ')
print('STATE 1: Binding Energies  [eV]')
for q in range(lenE):
    s = q+L1+1;   print('n = %0.0f' %s + '  EB = %0.3f' %EB[q])
print('State 1: n1 = %0.0f ' %n1 +' EB1 = %0.3f eV ' %EB1)
print('  ')

# STATE 2
HM2 = KM + UM2
# Eigenvalues [J] and eigenfunctions (eigenvectors)
ev, ef = eigsh(HM2, which="SM", k = num)
E = real(ev[ev<0]/se)                 # negative eigenvalues [eV]
lenE = len(E)
EB = -E
N2 = n2-L2-1
EB2 = EB[N2]

psi = zeros([N,len(E)]); psi2 = zeros([N,len(E)])
for c in range(len(E)):
    psi[:,c] = ef[:,c]
    psi2[:,c] = psi[:,c]**2
    area = simps(psi2[:,c],r)
    psi[:,c] = psi[:,c]/sqrt(area)
probD2 = psi**2    # probability density [1/m]
g2 = psi[:,N2]

for q in range(lenE):
    s = q+L2+1;   print('n = %0.0f' %s + '  EB = %0.3f' %EB[q])
print('State 2:  n2 = %0.0f ' %n2 +'  EB2 = %0.3f eV ' %EB2)

# Radial dipole moment contribution
fR = r*g1*g2 
IR = simps(fR,r)
print(' ')
print('Radial equations  IR = %0.3e ' %IR)
print('  ')


#%% ELECTRIC DIPOLE MOMENTS
px = e*abs(IR*IASx)
py = e*abs(IR*IASy)
pz = e*abs(IR*IASz)

p = sqrt(px**2 + py**2 + pz**2)

print(' ')
print('px = %0.3e ' %px + 'py = %0.3e ' %py  +   '   pz = %0.3e '%pz 
      + '   p = %0.3e' %p) 

#%% OSCILLATIONS   TRANSITION RATE   LIFETIME
f = abs(EB1-EB2)*se/h
dE = h*f/se
wL = cL/f
R = 16*pi**3*f**3*p**2/(3*eps0*h*cL**3)
tavg = 1/R
print(' ')
s = wL/sx; print('photon   f = %0.1e Hz' %f + '   lambda = %0.0f nm' %s
   + '  dE = %0.2f eV' %dE            )
print('Transition rate   R = %0.2e /s' %R)
s = tavg*1e9; print('Lifetime       tavg = %0.2f ns' %s)



#%% SAVE FIGURES
# fig1.savefig('a1.png')
fig2.savefig('a2.png')



#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)





