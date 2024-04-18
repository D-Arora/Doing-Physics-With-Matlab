# -*- coding: utf-8 -*-
"""
# qm004.py        April 2024
# Ian Cooper         matlabvisualphysics@gmail.com
# QUANTUM MECHANICS
#   EXPECTATION VALUES AND THE UNCERTAINTY PRINCIPLE
# # Website: https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm002.pdf
"""

import time
import math
import numpy as np
import pylab as py
from matplotlib import pyplot as plt 
from numpy import linspace,sin,cos,exp, zeros, ones, pi, diag, sqrt, real, imag
from scipy.integrate import odeint, quad, dblquad, simps

from matplotlib.ticker import (MultipleLocator, 
                               FormatStrFormatter, 
                               AutoMinorLocator) 

#%% FUNCTIONS

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


# INPUTS
# Quantum number
n = 2
# Grid points
N = 199
# Well width
a = 1e-9
# m to nm
L = 1e9

# SETUP
hbar = 1.054571817e-34; me = 9.11e-31; e = 1.602e-19; h = 6.626e-34
# X grid
x1 = 0; x2 = a; x = linspace(x1,x2,N); dx = x[2] - x[1]
# Wavefunction
y = sqrt(2/a)*sin(n*pi*x/a)

wL = 2*a/n
p = h/wL
K = p**2/(2*me)

# Probability
fn = y*y
prob = simps(fn,x)

# Expectation value & Uncertainty x
fn = y*x*y
xavg = simps(fn,x)
fn = y*x**2*y
x2avg = simps(fn,x)
deltax = sqrt(x2avg - xavg**2)      

# Expectation value & Uncertainty p
y1dash = firstDer(N,dx)@y
fn = y*y1dash
pavg = -1j*hbar*simps(fn,x)

y2dash = secondDer(N,dx)@y
fn = y*y2dash
p2avg = -hbar**2*simps(fn,x)
deltap = sqrt(p2avg - imag(pavg)**2)

delta = deltax*deltap

# Expectation KE
Kavg = -hbar**2*simps(fn,x)/(2*me)

# Console output
v = prob; print('prob =  %2.3f'  % v)
v = xavg*L; print('xavg =  %2.3f  nm'  % v)
v = deltax; print('deltax =  %2.3e  m'   % v)
v = deltap; print('deltap =  %2.3e  N.s'   % v)
v = delta; print('Delta =  %2.3e J.s'   % v)
v = Kavg/e; print('<K> =  %2.3e eV'   % v)
v = wL*L; print('lambda =  %2.3f nm'   % v)
v = p; print('p =  %2.3e N.s'   % v)
v = K/e; print('K =  %2.3e eV'   % v)

#%% GRAPHICS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)

fig, ax = plt.subplots(1)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('y  x10$^4$ ',color= 'black')
ax.set_xlabel('x [nm]',color = 'black')
fig.tight_layout()
ax.plot(x,y/1e4,'b', lw = 2)


# fig.savefig('a1.png')

