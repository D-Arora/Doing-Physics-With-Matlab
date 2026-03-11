# -*- coding: utf-8 -*-
'''

mnsET02.py.py      MARCH 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

CABLE EQUATION: Finite difference solution of cable equation

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mnsET01.pdf

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

#%% Model parameter
L = 1e-3
tau = 10e-3
C1 = 0.1
tMax = 12*tau
xMax = 5*L
xMin = -xMax
v0 = 0.1

nX = 199
x = linspace(xMin,xMax,nX)
dx = x[2] - x[1]

nC = int(nX/2) 
dt = C1*tau*dx**2/L**2
C2 = dt/tau

t = np.arange(0,tMax,dt)
nT = len(t)

s = 0.0001
V0 = v0*exp(-t**2/s)

v = zeros([nT,nX]) 
 
for T in range(nT-1):
    for X in range(2,nX-1):
         #v[T,nC] = v0
         v[T,nC] = V0[T] 
         v[T+1,X] = v[T,X] + C1*(v[T,X+1] - 2*v[T,X] + v[T,X-1]) - C2*v[T,X] 

#%%   FIG 1: Stimulus   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t [ms]'); ax.set_ylabel('v$_{m}$(t,0) [ mV ] ')
ax.grid()

ax.plot(t*1e3,V0*1e3,'b',lw = 2) 

fig1.tight_layout()
fig1.savefig('a1.png')

#%%   FIG 2
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)

def Graph2(q):
    yP = 1e3*v[q,:]; Q = t[q]/tau; Qstr = '%0.2f' %Q
    ax.plot(xP,yP,lw = 2, label = Qstr)
    
ax.set_xlabel('x / $\lambda$'); ax.set_ylabel('v$_m$ [ mV ] ')
ax.grid()
xP= x/L

A = array([100,300,500,700,1000,20e3,40e3],dtype=int)
for q in A:
   Graph2(q)

ax.legend(title = r't / $\tau$', loc = 'right', fontsize = 10)
ax.set_xlim([-2,2])
fig2.tight_layout()
fig2.savefig('a2.png')

#%%   FIG 3:    
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)

def Graph(q):
    yP = 1e3*v[R,q]; Q = x[q]/L; Qstr = '%0.2f' %Q
    ax.plot(xP,yP,lw = 2, label = Qstr)
    
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms]'); ax.set_ylabel('v$_m$ [ mV ] ')

ax.grid()

R = range(0,20000)
#R = range(0,nT-1)

xP = t[R]*1e3

A = [110,120,130,150,170]
for q in A:
   Graph(q)

ax.legend(title = 'x / $\lambda$', loc = 'right', fontsize = 10)
ax.set_xlim([0,30])
fig3.tight_layout()
fig3.savefig('a3.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


