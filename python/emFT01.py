# -*- coding: utf-8 -*-
"""
emPW01.py          mar 2025

COMPUTATIONAL OPTICS: SUPERPOSITION OF PLANE WAVES
     pLANES WAVES TRAVEL IN THE +z DIRECTION

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/emPW01.htm
"""

#%%
import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array, real, imag, conj 
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation, PillowWriter 
from scipy.integrate import simps

#%%   
# SETUP
nT = 999; nF = 999
f0 = 1e9; w0 = 2*pi*f0; T0 = 1/f0
t1 = -10*T0; t2 = 10*T0
f1 = 0*f0; f2 = 1.5*f0
E0 = 1
s = 4*T0
t = linspace(t1,t2,nT)
f = linspace(f1,f2,nF)
w = 2*pi*f

# Electric field and intensity
Et = E0*exp(-t**2/(2*s**2)) * exp(-1j*w0*t)
Et = Et + (0.8*E0)*exp(-t**2/(2*(2*s)**2)) * exp(-1j*w0*t/2)
SE = real(Et)**2
SEE = real(conj(Et)*Et)

# Fourier Transform
H  = zeros(nF) + 1j*zeros(nF)
for q in range(nF):
    h = Et * exp(1j*w[q]*t)
    HR = simps(real(h),t)
    HI = simps(imag(h),t)
    H[q]  = HR + 1j*HI 

Ef = H/np.sqrt(2*pi)
Sf = real(conj(Ef)*Ef)
S = Sf/max(Sf)

# Power
Pt = simps(SEE,t)    
Pf = simps(Sf,w)


#%%  [1D]    GRAPHICS 
plt.rcParams["figure.figsize"] = (7,3)
fig1, ax = plt.subplots(nrows=1, ncols=2)
q = s*1e9; fig1.suptitle('s = %0.0f  ns' % q)
ax[0].set_xlabel('t  [ns]',fontsize = 12)
ax[0].set_ylabel('E$_t$ [ a.u ] ',fontsize = 12)
ax[0].grid()
xP = t*1e9; yP = Et
ax[0].plot(xP,yP,'b',lw = 2)

ax[1].set_xlabel('t  [ns]',fontsize = 12)
ax[1].set_ylabel('S$_t$ [ a.u ] ',fontsize = 12)
ax[1].grid()
xP = t*1e9; yP = SE
ax[1].plot(xP,yP,'b',lw = 2)
yP = SEE; ax[1].plot(xP,yP,'r',lw = 2)
fig1.tight_layout()

plt.rcParams["figure.figsize"] = (7,3)
fig2, ax = plt.subplots(nrows=1, ncols=2)
q = s*1e9; fig2.suptitle('s = %0.0f  ns' % q)
ax[0].set_xlabel('f  [Hz]',fontsize = 12)
ax[0].set_ylabel('E$_f$ [ a.u ] ',fontsize = 12)
ax[0].grid()
xP = f; yP = H
ax[0].plot(xP,yP,'b',lw = 2)

ax[1].set_xlabel('f  [Hz]',fontsize = 12)
ax[1].set_ylabel('S$_f$ [ a.u ] ',fontsize = 12)
ax[1].grid()
xP = f; yP = Sf
ax[1].plot(xP,yP,'b',lw = 2)
fig2.tight_layout()


#%%  CONSOLE OUTPUT
q = Pt*1e9; print('Pt = %0.3f nW' %q)
q = Pf*1e9; print('Pf = %0.3f nW' % q)


#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')