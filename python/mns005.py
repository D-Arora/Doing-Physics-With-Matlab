# -*- coding: utf-8 -*-
'''
mns005.py      Jan 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

GATE VARIABLES m n h AND ION CHANNELS K+ Na+

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
# Resting value for membrane potential  [-65e-3  mV]
Vr = -65
# Membrane potential  [-100   50  mV]
Vmin = -100;   Vmax = 50;
# Reversal potentials [Na+  50    K+ -77   mV] 
VNa = 50;       
VK = -77;       
# Temperature  [20 degC]
T = 20
# Number of grid points  ( N such that dV = Vm - Vr <> -10e-3 ) 
N = 510
  
eps = 1e-16

#%%
# Membrane potential  [V]
Vm = linspace(Vmin,Vmax,N);
# Initialize arrays
An = zeros(N)
Bn  = zeros(N)
Am = zeros(N)
Bm  = zeros(N)
Ah = zeros(N)
Bh  = zeros(N)

V = 0

for c in range(N):
    [ An[c], Bn[c], Am[c], Bm[c], Ah[c], Bh[c] ] = AB(Vm[c])


#%% Gating variables m,h and n: asymptotic (steady-state) values
nINF = An/(An + Bn)
tauN = 1/(An + Bn)

mINF = Am/(Am + Bm)
tauM = 1/(Am + Bm)

hINF = Ah/(Ah + Bh)
tauH = 1/(Ah + Bh)


#%%    FIG 1:  Rate constants alpha (A) and beta(B) 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=2)

C = 0 
ax[C].set_xlabel('V$_m$  [ mV ]')
ax[C].set_ylabel(r'$\alpha$$_m$',fontsize = 14)
ax[C].set_title('sodium ')
ax[C].grid()
ax[C].set_xlim([-100,50])
ax[C].plot(Vm,Am,'b', lw = 2)

C = 1 
ax[C].set_xlabel('V$_n$  [ mV ]')
ax[C].set_ylabel(r'$\alpha$$_n$',fontsize = 14)
ax[C].set_title('potasium ')
ax[C].grid()
ax[C].set_xlim([-100,50])
ax[C].plot(Vm,An,'r', lw = 2)

fig1.tight_layout()
fig1.savefig('a1.png')

#%%    FIG 2:  Rate constants alpha (A) and beta(B) 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig2, ax = plt.subplots(nrows=1, ncols=2)

C = 0 
ax[C].set_xlabel('V$_m$  [ mV ]')
ax[C].set_ylabel(r'$\beta$$_n$',fontsize = 14)
ax[C].set_title('sodium ')
ax[C].grid()
ax[C].set_xlim([-100,50])
ax[C].plot(Vm,Bm,'b', lw = 2)
           
C = 1 
ax[C].set_xlabel('V$_n$  [ mV ]')
ax[C].set_ylabel(r'$\beta$$_n$',fontsize = 14)
ax[C].set_title('potasium ')
ax[C].grid()
ax[C].set_xlim([-100,50])
ax[C].plot(Vm,Bn,'r', lw = 2)

fig2.tight_layout()
fig2.savefig('a2.png')

#%%    FIG 3:  Rate constants alpha (A) and beta(B) 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=2)

C = 0 
ax[C].set_xlabel('V$_m$  [ mV ]')
ax[C].set_ylabel(r'$\alpha$$_h$',fontsize = 14)
ax[C].set_title('sodium ')
ax[C].grid()
ax[C].set_xlim([-100,50])
ax[C].plot(Vm,Ah,'b', lw = 2)
           
C = 1 
ax[C].set_xlabel('V$_n$  [ mV ]')
ax[C].set_ylabel(r'$\beta$$_h$',fontsize = 14)
ax[C].set_title('sodium ')
ax[C].grid()
ax[C].set_xlim([-100,50])
ax[C].plot(Vm,Bh,'b', lw = 2)

fig3.tight_layout()
fig3.savefig('a3.png')


#%%   FIG 4: Gating variables m,h and n: asymptotic (steady-state) values 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,7)
fig4, ax = plt.subplots(nrows=3, ncols=2)

R,C = 0, 0 
ax[R,C].set_xlabel('V$_m$  [ mV ]')
ax[R,C].set_ylabel(r'm$_{inf}$',fontsize = 14)
ax[R,C].set_title('sodium ')
ax[R,C].grid()
ax[R,C].set_xlim([-100,50])
ax[R,C].plot(Vm,nINF,'b', lw = 2)
           
R,C = 1, 0 
ax[R,C].set_xlabel('V$_m$  [ mV ]')
ax[R,C].set_ylabel(r'h$_{inf}$',fontsize = 14)
ax[R,C].set_title('sodium ')
ax[R,C].grid()
ax[R,C].set_xlim([-100,50])
ax[R,C].plot(Vm,hINF,'b', lw = 2)

R,C = 2, 0 
ax[R,C].set_xlabel('V$_m$  [ mV ]')
ax[R,C].set_ylabel(r'n$_{inf}$',fontsize = 14)
ax[R,C].set_title('potasium ')
ax[R,C].grid()
ax[R,C].set_xlim([-100,50])
ax[R,C].plot(Vm,nINF,'r', lw = 2)


R,C = 0, 1 
ax[R,C].set_xlabel('V$_m$  [ mV ]')
ax[R,C].set_ylabel(r'$\tau$$_n$',fontsize = 14)
ax[R,C].set_title('sodium ')
ax[R,C].grid()
ax[R,C].set_xlim([-100,50])
ax[R,C].plot(Vm,tauN,'b', lw = 2)
           
R,C = 1, 1 
ax[R,C].set_xlabel('V$_m$  [ mV ]')
ax[R,C].set_ylabel(r'$\tau$$_h$',fontsize = 14)
ax[R,C].set_title('sodium ')
ax[R,C].grid()
ax[R,C].set_xlim([-100,50])
ax[R,C].plot(Vm,tauH,'b', lw = 2)

R,C = 2, 1 
ax[R,C].set_xlabel('V$_m$  [ mV ]')
ax[R,C].set_ylabel(r'$\tau$$_n$',fontsize = 14)
ax[R,C].set_title('potasium ')
ax[R,C].grid()
ax[R,C].set_xlim([-100,50])
ax[R,C].plot(Vm,tauN,'r', lw = 2)

fig4.tight_layout()
fig4.savefig('a4.png')   



#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


