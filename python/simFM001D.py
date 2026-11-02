# -*- coding: utf-8 -*-
"""
Oct 2026         simFM001B.py

VIBRATIONS:
    FOURIER METHODS - FOURIER SERIES of a triangular wave

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/simFM001.pdfy


https://vcfw.org/pdf/Department/Physics/Fourier_series_python_code.pdf

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time
from scipy.integrate import solve_ivp, odeint, simpson
from numpy import pi, linspace, zeros, array, ones, sqrt, cos, sin, exp

tStart = time.time()
plt.close('all')

#%%
# INPUTS
# Number of time steps
Nt = 999
# Number of fourier components to sum
N = 31
# time span
# fundamental period
T1 = 1
# amplitude of wave
h = 1.5

#%% SETUP PARAMETERS
# Time mesh
t = linspace(0,T1,Nt)
# Frequency
f1 = 1/T1; w1 = 2*pi*f1

# width: sawtooth (1) or triangluar (0.5)
width = 0.5
#  phase angle shift in S
if width == 0.5: S = w1*t+pi/2
if width == 1.0: S = w1*t+pi
# Sawtooth or triangular function  fn
fn = h*signal.sawtooth(S,width)

#%% CALCULATE FOURIER COEFFICIENTS   n harmi=onic number
an = zeros(N); bn = zeros(N)
for n in range(1,N):
    fnC = fn * cos(n*w1*t); fnS = fn * sin(n*w1*t)
    an[n] = (2/T1)*simpson(fnC,t)
    bn[n] = (2/T1)*simpson(fnS,t)

#%% CALCULATE FOURIER SERIES fnFS
fnFS = zeros(Nt)
for c in range(1, N):
     fnFS = fnFS + an[c] * cos(c*w1*t) + bn[c] * sin(c*w1*t) 

#%% 1 WAVE
xP = t
txtX = 't [ s ]'; txtY = 'f (t)'; txt = ('n$_F$ = %0.0f' % N)
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.set_title(txt)
ax.grid()
yP = fn; ax.plot(xP,yP,'b',lw = 3,label = 'Orginal fn')
yP = fnFS; ax.plot(xP,yP,'r',lw = 2, label = 'Fourier series')
ax.axhline(0, color='black', linewidth=0.5)
ax.legend(fontsize = 10)
fig1.tight_layout()
fig1.savefig('a1.png')

#%% 2 Fourier coeffiicents
plt.rcParams['font.size'] = 12
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
xP = np.arange(1,12,1)
yP = bn[1:12]/bn[1]
markerline, stemlines, baseline =plt.stem(xP, yP, linefmt='b', markerfmt='bo')
plt.setp(markerline, markersize=7, markerfacecolor='b',label = 'b$_n$')
yP = an[1:12]
markerline, stemlines, baseline = plt.stem(xP, yP, linefmt='r', markerfmt='r*')
plt.setp(markerline, markersize=7, markerfacecolor='red',label = 'a$_n$')
ax.axhline(0, color='black', linewidth=1)
txtX = 'harmonic number n'; txtY = 'Fourier coeff. '
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.legend()
ax.grid()
fig2.tight_layout()
fig2.savefig('a2.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)

