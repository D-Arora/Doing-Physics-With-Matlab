# -*- coding: utf-8 -*-
'''
mns007.py      Jan 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

GEOMETRICAL ANAYLSIS OF [1D] DYNAMICAL SYSTEMS
   LEAK / FAST Na+ Channels

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mns007.pdf


'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.linalg import eig
import time

tStart = time.time()
plt.close('all')
    
N = 9999
#I0 = 0.6e-3
I0 = 0.88e-3

# Model parameters
CM = 10e-6      # Membrane capacitance [uF]
GL = 19e-3      # leak conductance [S]
GNa_max = 74e-3    # Na+ conductance  [S]
EL  = -67e-3     # Leak reversal potential / Nernst potential  [V]
ENa = 60e-3      # Na+ reversal potential
Vh = 19e-3       # m_inf
k  = 9e-3        # m_inf


#%% Na+ conductance as a function of V = VM   [V]
V1 = -100e-3; V2 = 50e-3; V = linspace(V1,V2,N)
m_inf = 1/(1 + np.exp((Vh - V)/k))
GNa = GNa_max * m_inf
# vDot  [V/s]

M = 399
v1 = zeros(M); v2 = zeros(M); v3 = zeros(M)
I1 = 0; I2 = 1.5e-3
Iext = linspace(I1,I2,M)

for m in range(M):
    Vdot = (Iext[m] - GL*(V - EL) - GNa*(V - ENa) ) / CM
    z = zeros(3); p = 0; 
    for c in range(N-1):
        q = Vdot[c]*Vdot[c+1]
        if q <= 0:
           z[p] = c
           p = int(p+1)     
    z = z.astype(int)
    Vzero = V[z]*1e3
    v1[m] = Vzero[0]
    v2[m] = Vzero[1]
    v3[m] = Vzero[2]


#%%   FIG 1: Bifurcation diagram 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel('V$_{SS}$  [ mV ]')
ax.set_xlabel('| I$_{ext}$ |  [ mA ]')
ax.grid()
ax.set_ylim([-80,60])
ax.set_xlim([0,1.5])
ax.plot(Iext*1e3,v1,'bo',ms= 1)
ax.plot(Iext*1e3,v2,'ro',ms= 1)
ax.plot(Iext*1e3,v3,'mo',ms= 1)
ax.plot([0.034,0.034],[-80,60],'g',lw = 1)
ax.plot([0.89,0.89],[-80,60],'g',lw = 1)
ax.text(0.2,-70,'stable',color = 'b')
ax.text(0.2,0,'unstable',color = 'r')
ax.text(0.2,40,'stable',color = 'm')
ax.text(1,30,'stable',color = 'b')
fig1.tight_layout()
fig1.savefig('a1.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


