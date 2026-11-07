# -*- coding: utf-8 -*-
"""
Nov 2026         simFM003.py

VIBRATIONS:
    FOURIER METHODS - FOURIER TRANSFORM - 
        SUPERPOSITION SINUSOIDAL FUNCTIONS

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/simFM003.pdf


"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time
from scipy.integrate import solve_ivp, odeint, simpson
from numpy import pi, linspace, zeros, array, ones, sqrt, cos, sin, exp, real, imag

tStart = time.time()
plt.close('all')

#%% TIME DOMAIN
Nt = 1999
t1 = 0; t2 = 5
t = linspace(t1,t2,Nt) 
A = array([1,1/3,1/5,1/7])
f0 = 5; f = array([1,3,5,7])
w = 2*pi*f

h = A[0]*sin(w[0]*t) + A[1]*sin(w[1]*t)+ A[2] *sin(w[2]*t) + A[3]*sin(w[3]*t)

#%%  FREQUENCY DOMAIN
Nf = 3999
f1 = -12; f2 = 12
f = linspace(f1,f2,Nf)
p  = 1j*2*pi
HR = zeros(Nf); HI = HR; H = HR+HI*1j
# Fourier transform
for c in range(Nf):
     g = h*exp(p*f[c]*t)
     H[c] = simpson(real(g),t) + 1j*simpson(imag(g),t)

psd = 2*abs(H)**2              # one-sided power spectral density                 

phase = np.angle(H)

powerF = simpson(abs(H)**2,f)
powerT = simpson(abs(h)**2,t)

#%% CONSOLE DISPLAY
print('Power (time domain) Pt =  %0.3f' % powerT)
print('Power (frequency domain) Pf =  %0.3f' % powerF)


#%% GRAPH DISPLAY
xP = t
txtX = 't [ s ]'; txtY = 'h(t)'
plt.rcParams['font.size'] = 10
fig1, ax = plt.subplots(nrows=3, ncols=1,figsize = (6,6))
C = 0
ax[C].set_xlabel(txtX)
ax[C].set_ylabel(txtY)
ax[C].grid()
yP = h; ax[C].plot(xP,yP,'b',lw = 2)

C = 1
xP = f
txtX = 'f [ Hz ]'; txtY = '| H(f) |'
ax[C].set_xlabel(txtX)
ax[C].set_ylabel(txtY)
ax[C].grid()
yP = abs(H); ax[C].plot(xP,yP,'b',lw = 2)

C = 2
xP = f
txtX = 'f [ Hz ]'; txtY = 'psd(f)'
ax[C].set_xlabel(txtX)
ax[C].set_ylabel(txtY)
ax[C].grid()
xP = f[f>0]; yP = psd[f>0]; ax[C].plot(xP,yP,'b',lw = 2)

fig1.tight_layout()
fig1.savefig('a1.png')

#%% Phase
txtX = 'f [ Hz ]'; txtY = 'phase  [ rad ]'
plt.rcParams['font.size'] = 10
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
xP = f; yP = phase; ax.plot(xP,yP,'b',lw = 1)

fig2.tight_layout()
fig2.savefig('a2.png')

#%%   Fourier transform  real and imaginary parts
txtX = 'f [ Hz ]'; txtY = 'H(f) '
plt.rcParams['font.size'] = 10
fig3, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
xP = f; yP = real(H); ax.plot(xP,yP,'b',lw = 1,label = 'real')
xP = f; yP = imag(H); ax.plot(xP,yP,'r',lw = 1, label = 'imag')
ax.legend()
fig3.tight_layout()
fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)


