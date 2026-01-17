# -*- coding: utf-8 -*-
'''
mns005.py      Jan 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

GATE VARIABLES m n h AND ION CHANNELS K+ Na+
Voltage clamp simulations

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/mns005.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt, exp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.linalg import eig
import time

tStart = time.time()
plt.close('all')


#%%
def AB(V):
    dV = (V - Vr)
    phi = 3**((T-6.3)/10)

    An = phi * (eps + 0.10 - 0.01 * dV) / (eps + exp(1 - 0.1 * dV) - 1)
    Am = phi * (eps + 2.5 - 0.1 * dV) / (eps + exp(2.5 - 0.1 * dV) - 1)
    Ah = phi * 0.07 * exp(-dV / 20)

    Bn = phi * 0.125 * exp(-dV / 80)
    Bm = phi * 4 * exp(-dV/18)
    Bh = phi * 1 / (exp(3.0 - 0.1 * dV) + 1)  
    return [An, Bn, Am, Bm, Ah, Bh]    

#%%
# Resting value for membrane potential  [mV]
Vr = -65
# Reversal potentials [mV] 
ENa = 50       
EK = -77   
EL = - 75.6    
# Temperature  [degC]
T = 6.3

gK_max   = 36         # max conductance K         [mS.cm-2]
gNa_max  = 120        # max conductance Na        [mS.cm-1]
gL_max   = 0.3        # max conductance leakage   [mS.cm-2]

# Number of grid points  ( N such that dV = Vm - Vr <> -10e-3 ) 
N = 999
eps = 1e-16
# Time span
t = linspace(0,20,N)
dt = t[2] - t[1]
# Clamp potential [mV]
dV = 80
V = Vr *ones(N)
V[t>2.5] = Vr + dV
V[t>12.5] = Vr

#%% # Initialize arrays
JNa = zeros(N); JK  = zeros(N)
gNa = zeros(N); gK = zeros(N)
m = zeros(N); h  = zeros(N); n = zeros(N)

An = zeros(N); Bn  = zeros(N)
Am = zeros(N); Bm  = zeros(N)
Ah = zeros(N); Bh  = zeros(N)

#%%

for c in range(N):
    [ An[c], Bn[c], Am[c], Bm[c], Ah[c], Bh[c] ] = AB(V[c])
    
n[0] = An[0] / (An[0] + Bn[0])
m[0] = Am[0] / (Am[0]+ Bm[0])
h[0] = Ah[0] / (Ah[0] + Bh[0])

for c in range(N-1):
    n[c+1] = n[c] + dt * ( An[c] * (1-n[c]) - Bn[c] * n[c] ) 
    m[c+1] = m[c] + dt * ( Am[c] * (1-m[c]) - Bm[c] * m[c] ) 
    h[c+1] = h[c] + dt * ( Ah[c] * (1-h[c]) - Bh[c] * h[c] ) 

#%% Conductances  mS.cm-2   A.mV-1.cm-2
gK  = gK_max  * n**4
gNa = gNa_max * m**3 * h

# Current densities  A --> mA
JK  = gK   * (V - EK) / 1000
JNa = gNa  * (V - ENa) /1000
JL = gL_max * (V - EL) / 1000
Jm = JK + JNa + JL


#%%    FIG 1:  VOLTAGE CLAMP
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms ]')
ax.set_title('Voltage clamp  dV = %0.0f mV '%dV +
             '  T = %0.1f $^o$C' % T, fontsize = 10)
ax.grid()
#ax.set_ylim([0,1.1])
ax.set_xlim([0,20])
ax.plot(t,V,'k',lw = 2)
fig1.tight_layout()
fig1.savefig('a1.png')

#%%    FIG 2:  GATING VARIABLES n m h
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms ]')
ax.set_title('Gating variables n m h',fontsize = 10)
ax.grid()
ax.set_ylim([0,1.1]); ax.set_xlim([0,20])
ax.plot(t,n,'r',lw = 2,label = 'n')
ax.plot(t,m,'b',lw = 2,label = 'm')
ax.plot(t,h,'m',lw = 2,label = 'h')
ax.legend(fontsize = 10)
fig2.tight_layout()
fig2.savefig('a2.png')

#%%    FIG 3:  GATING VARIABLES n m h
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms ]')
ax.set_ylabel('g  [ mS.cm$^{-2}$]')
ax.set_title('Conductance g',fontsize = 10)
ax.grid()
ax.set_ylim([0,40])
ax.set_xlim([0,20])
ax.plot(t,gK,'r',lw = 2,label = 'K$^+$')
ax.plot(t,gNa,'b',lw = 2,label = 'Na$^+$')
ax.legend(fontsize = 10)
fig3.tight_layout()
fig3.savefig('a3.png')

#%%    FIG 4:  CUUURENT DENSITIES
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms ]')
ax.set_ylabel('J  [ mA.cm$^{-2}$]')
ax.set_title('current densities J',fontsize = 10)
ax.grid()
#ax.set_ylim([0,40])
ax.set_xlim([0,20])
ax.plot(t,JK,'r',lw = 1,label = 'K$^+$')
ax.plot(t,JNa,'b',lw = 1,label = 'Na$^+$')
ax.plot(t,JL,'m',lw = 1,label = 'J$_L$')
ax.plot(t,Jm,'k',lw = 2,label = 'J$_m$')
ax.legend(fontsize = 10, ncols = 4, loc = 'lower right',frameon=False)
fig4.tight_layout()
fig4.savefig('a4.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


