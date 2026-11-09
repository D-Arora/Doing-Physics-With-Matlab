# -*- coding: utf-8 -*-

"""
Nov 2026         simFM004A.py

FOURIER METHODS - FOURIER TRANSFORM  
  ECG RECORDING

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/simFM004.pdf

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time
from scipy.integrate import solve_ivp, odeint, simpson
from numpy import pi, linspace, zeros, array, ones, sqrt, cos, sin, exp, real, imag
import scipy.io

from scipy.io import wavfile

tStart = time.time()
plt.close('all')

#%% TIME DOMAIN: Load the .mat file into a Python dictionary
ecg_data = scipy.io.loadmat('ecg.mat')
# Extract the ecg data
x = ecg_data['ecg']
# Time span [s]
t = linspace(0,18,len(x))
# ecg data function
h = x[:,0]

#%%  FREQUENCY DOMAIN
Nf = 3999
f1 = -19; f2 = 19
f = linspace(f1,f2,Nf)
p  = 1j*2*pi
HR = zeros(Nf); HI = HR; H = HR+HI*1j
# Fourier transform
for c in range(Nf):
     g = h*exp(p*f[c]*t)
     H[c] = simpson(real(g),t) + 1j*simpson(imag(g),t)
# one-sided power spectral density 
psd = 2*abs(H)**2                   

#%% PLOT 1: ecg recording
txtX = 't [ s ]'; txtY = 'ECG) '
plt.rcParams['font.size'] = 10
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (6,3))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
xP = t; yP = h; ax.plot(xP,yP,'b',lw = 1)
ax.set_xticks(range(0,19,2))
fig1.tight_layout()
fig1.savefig('a1.png')

#%% PLOT 2: Fourier transform   H
txtX = 'f [ Hz ]'; txtY = '| H(f) '
plt.rcParams['font.size'] = 10
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (6,3))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
#xP = f; yP = abs(H); ax.plot(xP,yP,'b',lw = 1)
R = np.arange(2050,3990,1) 
xP = f[R]; yP = abs(H[R]); ax.plot(xP,yP,'b',lw = 1)
ax.set_xticks(range(0,19,1))
ax.set_xlim([0,18])
fig2.tight_layout()
fig2.savefig('a2.png')

#%% PLOT 2: Fourier transform   psd
txtX = 'f [ Hz ]'; txtY = 'psd '
plt.rcParams['font.size'] = 10
fig3, ax = plt.subplots(nrows=1, ncols=1,figsize = (6,3))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
#xP = f; yP =psd; ax.plot(xP,yP,'b',lw = 1)
R = np.arange(2050,3990,1) 
xP = f[R]; yP = psd[R]/max(psd[R]); ax.plot(xP,yP,'b',lw = 1)
ax.set_xticks(range(0,13,1))
ax.set_xlim(0,12)
fig3.tight_layout()
fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)
