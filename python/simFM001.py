# -*- coding: utf-8 -*-
"""
Oct 2026         simFM001.py

VIBRATIONS:
    FOURIER METHODS - FOURIER SERIES of a square wave

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/simFM001.pdfy

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, odeint, simpson
from numpy import pi, arange, linspace, zeros, array, ones, sqrt, exp, cos, sin 
import time

tStart = time.time()
plt.close('all')

#%%
# Frequency of square wave pulse [Hz]
f0 = 0.5
# Number of grid points
N = 999
# amplitude of square wave pulse [a.u.]
h = 1.2
# Time span for square wave / time mesh  [s]
t1 = -4; t2 = 12
t = linspace(t1,t2,N)
w0 = 2*pi*f0
# Sine function / square wave pulse
h0 = 0.750      # DC off-set
fSine = sin(w0*t) + h0
fSquare = fSine
fSquare[fSine > h0] = h+h0
fSquare[fSine < h0] = -h+h0


#%% Fourier series
n = 31
fn = zeros(N)
for c in range(n):
    co = 2*c+1
    fn = fn + (1/co)*sin(w0*co*t)
       
fn = h0 + (4*h/pi)*fn
    
#%% FOURIOR ANALYSIS
p  = 1j*2*pi
Fmax = 6
Fmin = 0.01
nF = 999
F = np.linspace(Fmin,Fmax,nF)
HR = np.zeros(nF); HI = HR; H = HR; H = HR+HR*1j
for c in range(nF):
     g = fSquare*np.exp(p*F[c]*t)
     H[c] = simpson(g,t)


psd = 2*np.conj(H)*H   
psd = np.real(psd/max(psd))   

Hn = sqrt(psd)

Find = np.where(np.real(psd)>= 0.98)   #0.98
II = Find[-1]
II = II[-1]
Fpeak = abs(F[II])

#%% 1 Square pulse
xP = t; yP = fSquare
txtX = 't [ s ]'; txtY = 'f (t)'
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
ax.plot(xP,yP,'b',lw = 2)
fig1.tight_layout()
fig1.savefig('a1.png')

#%% 2 Fourier Series
xP = t; yP = fSquare
txtX = 't [ s ]'; txtY = 'f (t)'
txt = ('n = %0.0f' %n)
plt.rcParams['font.size'] = 10
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.set_title(txt)
ax.grid()
yP = fn;      ax.plot(xP,yP,'r',lw = 2)
yP = fSquare; ax.plot(xP,yP,'b',lw = 1)
fig2.tight_layout()
fig2.savefig('a2.png')

#%% 3 Fourier transform
plt.rcParams['font.size'] = 12
fig3, ax = plt.subplots(nrows=3, ncols=1,figsize = (5,7))

xP = F/Fpeak; yP = Hn
txtX = r'f / f$_{peak}$'; txtY = '| H | '
ax[0].plot(xP,yP,'b',lw = 2)
ax[0].set_xlabel(txtX)
ax[0].set_ylabel(txtY)
ax[0].set_xticks(array([1,3,5,7,9,11]))
ax[0].grid()

xP = F/Fpeak; yP = psd
txtX = r'f / f$_{peak}$'; txtY = 'psd '
ax[1].plot(xP,yP,'b',lw = 2)
ax[1].set_xlabel(txtX)
ax[1].set_ylabel(txtY)
ax[1].set_xticks(array([1,3,5,7,9,11]))
ax[1].grid()

xP = F/Fpeak; yP = 10*np.log10(psd)
txtX = r'f / f$_{peak}$'; txtY = 'psd  [ dB ]'
ax[2].plot(xP,yP,'b',lw = 2)
ax[2].set_xlabel(txtX)
ax[2].set_ylabel(txtY)
ax[2].set_xticks(array([1,3,5,7,9,11]))
ax[2].grid()

fig3.tight_layout()
fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)

