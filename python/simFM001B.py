# -*- coding: utf-8 -*-
"""
Oct 2026         simFM001B.py

VIBRATIONS:
    FOURIER METHODS - FOURIER SERIES of a sawtooth wave

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
N = 51
# time span
t1 = 0; t2 = 2
# fundamental period
T1 = 0.5
# amplitude of wave
h = 1.5

# SETUP PARAMETERS
# Time mesh
t = linspace(t1,t2,Nt)
# Frequency
f1 = 1/T1; w1 = 2*pi*f1
# sawtooth phase angle S
S = w1*t+pi
# width: sawtooth (1) or triangluar (0.5)
width = 1
# Sawtooth function  fn
fn = h*signal.sawtooth(S,width)

#%% CALCULATE FOURIER SERIES fnFS
fnFS = zeros(Nt)
for n in range(1, N+1):
    # Calculate the bn coefficient: 2 * (-1)**(n+1) / n
    bn = (2/pi) * ((-1) ** (n + 1)) / n
    # Add the nth harmonic component to the sum
    fnFS += bn * sin(n*w1*t) 
fnFS = h*fnFS

#%% FOURIOR ANALYSIS
p  = 1j*2*pi
Fmax = 16
Fmin = 0.01
nF = 1999
F = np.linspace(Fmin,Fmax,nF)
HR = np.zeros(nF); HI = HR; H = HR; H = HR+HR*1j
for c in range(nF):
     g = fn*exp(p*F[c]*t)
     H[c] = simpson(g,t)

psd = 2*np.conj(H)*H   
psd = np.real(psd/max(psd))   

Hn = sqrt(psd)

Find = np.where(np.real(psd)>= 0.98)   #0.98
II = Find[-1]
II = II[-1]
Fpeak = abs(F[II])

#%% 1 WAVE
xP = t
txtX = 't [ s ]'; txtY = 'f (t)'; txt = ('n$_F$ = %0.0f' % N)
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.set_title(txt)
ax.grid()
yP = fn; ax.plot(xP,yP,'b',lw = 3)
yP = fnFS; ax.plot(xP,yP,'r',lw = 2)
ax.axhline(0, color='black', linewidth=0.5)
fig1.tight_layout()
fig1.savefig('a1.png')

#%% 3 Fourier transform
plt.rcParams['font.size'] = 12
fig3, ax = plt.subplots(nrows=3, ncols=1,figsize = (6,7))

xP = F/Fpeak; yP = Hn
txtX = r'f / f$_{peak}$'; txtY = '| H | '
ax[0].plot(xP,yP,'b',lw = 2)
ax[0].set_xlabel(txtX)
ax[0].set_ylabel(txtY)
#ax[0].set_xticks(array([1,3,5,7,9,11]))
ax[0].grid()

xP = F/Fpeak; yP = psd
txtX = r'f / f$_{peak}$'; txtY = 'psd '
ax[1].plot(xP,yP,'b',lw = 2)
ax[1].set_xlabel(txtX)
ax[1].set_ylabel(txtY)
#ax[1].set_xticks(array([1,3,5,7,9,11]))
ax[1].grid()

xP = F/Fpeak; yP = 10*np.log10(psd)
txtX = r'f / f$_{peak}$'; txtY = 'psd  [ dB ]'
ax[2].plot(xP,yP,'b',lw = 2)
ax[2].set_xlabel(txtX)
ax[2].set_ylabel(txtY)
#ax[2].set_xticks(array([1,3,5,7,9,11]))
ax[2].grid()

fig3.tight_layout()
fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)

