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
from scipy.io.wavfile import write
import os
import random

tStart = time.time()
plt.close('all')

#%% TIME DOMAIN
Nt = 9999
 
f0 = 500  
fnoise = 50
#Anoise = np.random.uniform(-1.0, 1.0, k = 5)
Anoise = zeros(Nt)
Anoise = [np.random.uniform(0.8, 1) for _ in range(Nt)]
t1 = 0; t2 = 20/f0
t = linspace(t1,t2,Nt)
h = sin(2*pi*f0*t) + Anoise*sin(2*pi*fnoise*t)

#%%  FREQUENCY DOMAIN
Nf = 3999
f1 = -600; f2 = 600
f = linspace(f1,f2,Nf)
p  = 1j*2*pi
HR = zeros(Nf); HI = HR; H = HR+HI*1j
# Fourier transform
for c in range(Nf):
     g = h*exp(p*f[c]*t)
     H[c] = simpson(real(g),t) + 1j*simpson(imag(g),t)

H[abs(f)<100] = 0
psd = 2*abs(H)**2              # one-sided power spectral density      


#%% INVERSE FOURIER TRANSFORM
Hinv = zeros(Nt)+1j*zeros(Nt)

#HRinv = zeros(Nf); HIinv = HRinv; Hinv = HRinv + HIinv*1j
for c in range(Nt):
     g = H*exp(-p*f*t[c])
     Hinv[c] = simpson(real(g),f) + 1j*simpson(imag(g),f)

#%% 1  GRAPH DISPLAY
txtX = 't [ ms ]'; txtY = 'h(t)'
plt.rcParams['font.size'] = 10
fig1, ax = plt.subplots(nrows=3, ncols=1,figsize = (6,6))
C = 0
ax[C].set_xlabel(txtX)
ax[C].set_ylabel(txtY)
ax[C].grid()
xP = t*1e3; yP = h; ax[C].plot(xP,yP,'b',lw = 2)

C = 1
xP = f
txtX = 'f [ Hz ]'; txtY = '| H(f) |'
ax[C].set_xlabel(txtX)
ax[C].set_ylabel(txtY)
ax[C].grid()
yP = abs(H)/max(abs(H)); ax[C].plot(xP,yP,'b',lw = 2)

C = 2
xP = f
txtX = 'f [ Hz ]'; txtY = 'psd(f)'
ax[C].set_xlabel(txtX)
ax[C].set_ylabel(txtY)
ax[C].grid()
xP = f[f>0]; yP = psd[f>0]/max(psd); ax[C].plot(xP,yP,'b',lw = 2)

fig1.tight_layout()
fig1.savefig('a1.png')


#%%  3 Fourier transform  real and imaginary parts
# txtX = 'f [ Hz ]'; txtY = 'psd  [dB] '
# plt.rcParams['font.size'] = 10
# fig3, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
# ax.set_xlabel(txtX)
# ax.set_ylabel(txtY)
# ax.grid()
# ax.set_ylim([-20,0])
# xP = f
# yP = 10*np.log10(psd/max(psd))
# ax.plot(xP,yP,'b',lw = 2)
# fig3.tight_layout()
# fig3.savefig('a3.png')

#%%   4 INVERSE FOURIER TRANSFORM
txtX = 't [ ms ]'; txtY = 'signal '
plt.rcParams['font.size'] = 10
fig4, ax = plt.subplots(nrows=1, ncols=1,figsize = (6,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
xP = t*1e3; yP = Hinv/max(Hinv); ax.plot(xP,yP,'b',lw = 2,label = 'inv')
xP = t*1e3; yP = h/max(h); ax.plot(xP,yP,'r',lw = 1, label = 'h')
ax.legend()
fig4.tight_layout()
fig4.savefig('a4.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)



