# -*- coding: utf-8 -*-
'''

mnsET01.py      MARCH 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

The Passive Cable Equation

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mnsIZH01.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, log, linspace, zeros, array, ones, sqrt, exp 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.linalg import eig
import time
import matplotlib.animation as animation

tStart = time.time()
plt.close('all')

#%%  STEADY-STATE SOLUTION OF CABLE EQUATION 

rm = 1.0      # membrane specific resistance  [ohm.m**2]
rL = 1.0      # specific longitudial resistance [ohm.m]
v0 = 0.1      # Initial steadt-state membrane potential [V]
a = array([2e-6,8e-6])     # cable radius [m]

L = (a*rm/rL)**0.5                 # lambda  [m]
Rinp = (rm*rL/2)/(pi*a**(3/2))     # [ohm]
Lh = -L*log(0.5)                   # half-space distance
Ginp = 1/Rinp

N = 599; xMax = 15e-3
x = linspace(0,xMax,N)

v1 = v0*exp(-x/L[0])
v2 = v0*exp(-x/L[1])
  
#%%   FIG 1: t vs V   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('x [ mm ]'); ax.set_ylabel('v$_{SS}$  [ mV ]')
ax.grid()
ax.set_xlim([0,15])
xP = x*1e3; yP = v1*1e3

q1 = a[0]*1e6; q2 = L[0]*1e3; q3 = Lh[0]*1e3; q4 = Ginp[0]*1e9
qstr = 'a = %0.2f $\mu$m' %q1 + \
 ' $\lambda$ = %0.2f mm' %q2 + \
 ' $\lambda$$_h$ = %0.2f mm' %q3 + \
 ' $G_{inp}$ = %0.2f nS' %q4

ax.plot(xP,yP,'b', lw = 2,label = qstr ) 
yP = v2*1e3; q = a[1]*1e6

q1 = a[1]*1e6; q2 = L[1]*1e3; q3 = Lh[1]*1e3; q4 = Ginp[1]*1e9
qstr = 'a = %0.2f $\mu$m' %q1 + \
 ' $\lambda$ = %0.2f mm' %q2 + \
 ' $\lambda$$_h$ = %0.2f mm' %q3 + \
 ' $G_{inp}$ = %0.2f nS' %q4


ax.plot(xP, yP,'r', lw = 2, label = qstr)
ax.plot([0,15],[50,50],'m',lw=1)
ax.legend(fontsize = 10)
fig1.tight_layout()
fig1.savefig('a1.png')


#%%  LINEAR MEMBRANE POTENTIAL PROFILE
V0 = 100
nX, xLim = 299, 10
x = linspace(-xLim,xLim,nX)
tm = 5
nT = 501
t = linspace(0,8*tm,nT)
g = V0*exp(-t/tm)


#%%   FIG 2:  
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t'); ax.set_ylabel('v$_m$ ')
ax.grid()
ax.plot(t,x[-1]*g,'r',lw = 2)
fig2.tight_layout()
fig2.savefig('a2.png')

#   FIG 3:    
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.grid()
T = array([0,20,50,100,150,500])
for c in range(6):
    vm = x*g[T[c]]
    ax.plot(x,vm, lw = 2,label = int(t[T[c]]))
ax.set_xlabel('x'); ax.set_ylabel('v$_m$ ')
ax.legend(fontsize = 10)
fig3.tight_layout()
fig3.savefig('a3.png')

         
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


