# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 18:57:07 2024

@author: Owner
"""

# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.signal import chirp, find_peaks, peak_widths
import time
from numpy.linalg import eig



sF = 1
q = np.pi
t1 = 0
t2 = 5
Fmin = -2
Fmax = 2
b = q/6
t = np.linspace(t1,t2,999)
x = np.sin(2*np.pi*sF*t+ b)
#x = np.exp(-np.pi*t**2)


#%% FOURIER TRANSFORM - frequency spectrum

def simpson1d(f,xMin,xMax):
    N = len(f)
    h = (xMax - xMin) / (N - 1)
    
    integral = (h/3) * (f[0] + 2*sum(f[:N-2:2]) \
            + 4*sum(f[1:N-1:2]) + f[N-1])
 
    if N%2 == 0:
        integral = 'N must be an odd number'
        print('integral')
    return integral

p  = 1j*2*np.pi

nF = 299
F = np.linspace(Fmin,Fmax,nF)
HR = np.zeros(nF); HI = HR; H = HR; HH = HR+HR*1j
for c in range(nF):
     g = x*np.exp(p*F[c]*t)
     gR = np.real(x*np.exp(p*F[c]*t))
     gI = np.imag(x*np.exp(p*F[c]*t))
     HR[c] = simpson1d(gR,t1,t2)
     HI[c] = simpson1d(gI,t1,t2)
     HH[c] = simpson1d(g,t1,t2)

#H = HR + HI*1j
psd = HR**2 + HI**2
#psd = 2*np.conj(H)*H   
psd = psd/max(psd)   

fig = plt.figure(figsize = (4, 3))
plt.plot(t,x)

fig = plt.figure(figsize = (4, 3))
#plt.plot(F,psd)
plt.plot(F,np.conj(HH)*HH)
plt.grid('visible')
