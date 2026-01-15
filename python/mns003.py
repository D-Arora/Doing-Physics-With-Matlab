# -*- coding: utf-8 -*-
'''
mns003.py      Jan 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

THE NEURON AS A CAPACITOR

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/mns003.pdf


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


#%%
def VDOT(t,V, J):
    fn = (-g*(V - E) + J)/c
    return fn    
    
#%% Time span
N = 999
tS = 0.2
t = linspace(0,tS,N)
h = t[2] - t[1]
# Model parameters
c = 1
E = -77
g = 36
Vreset = -100
Vrest = E
vTH = -50
# Memebrane potential:  V0 initial condition
V = zeros(N) 

V0 = -100

V[0] = V0

# External current Jext  
# No spike flag = 0   Spike flag = 1 
flag = 1
t1,t2 = 0.05, 0.10     # pulse duration dt = t2 - t1
Jmax = 1200            # pulse height
J = zeros(N)

if flag == 1:
   J[t>t1] = Jmax
   J[t>t2] = 0
   V0 = E
   V[0] = V0

for n in range(N-1):
    if flag == 1:
       if V[n-1] > vTH:
          V[n-1] = 0
          V[n] = Vreset
          V[0] = E
    k1 = VDOT(t[n],V[n],J[n])
    k2 = VDOT(t[n]+0.5*h,V[n]+k1*h/2,J[n])
    k3 = VDOT(t[n]+0.5*h,V[n]+k2*h/2,J[n])
    k4 = VDOT(t[n]+1.0*h,V[n]+k3*h,J[n]);
    V[n+1] = V[n] + h*(k1+2*k2+2*k3+k4)/6


Jm = g*(V - E)      # uA/cm^2
  
   
#%%   FIG 1: t vs V   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('V$_M$  [ mV ]')
ax.grid()
ax.set_ylim([-100,0])
ax.plot(t,V,'b', lw = 2)


if flag == 0:
  tau = c/g
  if V0 < Vrest:
     vC = V0 + (Vrest - V0)*(1-np.exp(-1))
     ax.set_ylim([-100,-70])
     ax.plot([tau,tau],[V0,vC],'r',lw=1)
     ax.plot([0,tau],[vC,vC],'r',lw=1)
     ax.text(0.045,-86,'tau = %0.3f ms' %tau)
     ax.text(0.045,-90,'v(tau) = %0.0f mV' %vC)
   
  if V0 > Vrest:
      ax.plot([0,tS],[Vrest,Vrest],'m',lw=1)
      vC = Vrest - (Vrest + V0)*np.exp(-1)
      #ax.set_ylim([0,-80])
      ax.plot([tau,tau],[Vrest,vC],'r',lw=1)
      ax.plot([0,tau],[vC,vC],'r',lw=1)
      ax.text(0.045,-40,'tau = %0.3f ms' %tau)
      ax.text(0.045,-50,'v(tau) = %0.0f mV' %vC)   

fig1.tight_layout()
# fig1.savefig('a1.png')

#%%   FIG 2: t vs J   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel(r'J$_m$  [ $\mu$A.cm$^{-2}$ ]')
ax.grid()
ax.plot(t,Jm,'b',lw = 2,label = 'J$_m$')
ax.plot(t,J,'r',lw = 0.5,label = 'J$_{ext}$')
ax.legend()
fig2.tight_layout()
# fig2.savefig('a2.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


