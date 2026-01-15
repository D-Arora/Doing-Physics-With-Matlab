# -*- coding: utf-8 -*-
'''
mns004.py      Jan 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

SPIKING NEURONS
    LEAKY INTEGRATE-AND-FIRE MODEL


Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/mns004.pdf


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
def VDOT(t,V, Iext):
    fn = (-(V-vRest) + R*Iext)/tau
    return fn    

#%% INPUTS AND MODEL PARAMETERS
# Initial membrane potential >>>
V0 = -60
# External current Iext >>>
flag = 3

flagP = 0      # flag for ramp input

# Time span
N = 9999
tS = 100
t = linspace(0,tS,N)
h = t[2] - t[1]

# Model parameters
R = 1e7
tau = 10
vReset = -80
vRest = -60
vTH = -30
vSpike = 0
# Membrane potential:  V0 initial condition
V = zeros(N) 
V[0] = V0
# Absolute refractrory period
tARP = 5     
nARP = round(tARP / h)
d = -1

#%%  EXTERNAL CURRENT STIMULUS 
Iext = zeros(N)
match flag:
    case 0:           # Zero external current
      Iext = zeros(N)
    case 1:           # Single  pulse
      Imax = 4.75e-6
      t1,t2 = 10, 20     # pulse duration dt = t2 - t1
      Iext[t>t1] = Imax
      Iext[t>t2] = 0

    case 2:           # Two pulses
      Imax = 4.8e-6
      t1,t2 = 10, 20     # pulse duration dt = t2 - t1
      Iext[t>t1] = Imax
      Iext[t>t2] = 0
      Imax = 7e-6
      t1,t2 = 22, 32     # pulse duration dt = t2 - t1
      Iext[t>t1] = Imax
      Iext[t>t2] = 0  
   
    case 3:           # Multiple pulses
        Imax = 4.9e-6
        T = 10
        y = sin(2*pi*t/T)
        Iext[y >= 0]: Iext = Imax
    
    case 4:          # Step input
        Imax = 10.0e-6
        T = 10
        Iext[t>T] = Imax
    
    case 5:         # Ramp input
        Iext = 20e-8*t   
        Imax = max(Iext)
        flagP = 1
        
#%% Runge-Kutta
for n in range(N-1):
    k1 = VDOT(t[n],V[n],Iext[n])
    k2 = VDOT(t[n]+0.5*h,V[n]+k1*h/2,Iext[n])
    k3 = VDOT(t[n]+0.5*h,V[n]+k2*h/2,Iext[n])
    k4 = VDOT(t[n]+1.0*h,V[n]+k3*h,Iext[n]);
    V[n+1] = V[n] + h*(k1+2*k2+2*k3+k4)/6
    
    if V[n+1] > vTH:
       V[n] = vSpike
       V[n+1] = vReset
       d = nARP
   
    if d > 0:
      V[n+1] = vReset
      d = d - 1

#%%   FIG 1: t vs V   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('V$_M$  [ mV ]')
ax.grid()
ax.set_ylim([-100,0])
VRest = ones(N)*vRest
ax.plot(t,VRest,'r', lw = 1)
ax.plot(t,vTH*ones(N),'r', lw = 1)

ax.plot(t,V,'k', lw = 2)

fig1.tight_layout()
fig1.savefig('a1.png')

#%%   FIG 2: t vs IEXT   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2)
fig2, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel(r'I  [ nA ]')
ax.grid()

ax.plot(t,Iext*1e6,'r',lw = 0.5,label = 'I$_{ext}$')

fig2.tight_layout()
fig2.savefig('a2.png')

#%%    FIG 3:  t vs V  and t Vs Iext 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=2, ncols=1)
 
ax[0].set_xlabel('t  [ ms ]'); ax[0].set_ylabel('V$_M$  [ mV ]')
q = 1e6*Imax; ax[0].set_title('I$_{max}$ = %0.0f nA'%q)
ax[0].grid()
ax[0].set_ylim([-100,0])

ax[0].plot(t,VRest*ones(N),'r', lw = 1)
ax[0].plot(t,vTH*ones(N),'r', lw = 1)
ax[0].plot(t,V,'b', lw = 2)

ax[1].grid()
ax[1].set_xlabel('t  [ ms ]'); ax[1].set_ylabel('I  [ nA ]')
ax[1].plot(t,Iext*1e6,'r',lw = 2,label = 'I$_{ext}$')

fig3.tight_layout()
fig3.savefig('a3.png')


#%%  FIG 4: ISI and firing rate as a function Iext

if flagP ==1:
   peaks, _ = find_peaks(V)
   tPeaks = t[peaks]
   L = len(tPeaks)
   ISI = (tPeaks[1:L] - tPeaks[0:L-1])*1e-3
   fp = 1/ISI
   Ip = Iext[peaks]

   plt.rcParams['font.size'] = 12
   plt.rcParams["figure.figsize"] = (5,2)
   fig4, ax = plt.subplots(nrows=1, ncols=1)
    
   ax.set_xlabel('I$_{ext}$  [ nA ]'); ax.set_ylabel('f  [ Hz]')
   ax.grid()
   ax.set_xlim([0,40])
   ax.plot(Ip[1:L]*1e6,fp,'b',lw = 2)

   fig4.tight_layout()
   fig4.savefig('a4.png')   
         
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


