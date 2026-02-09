# -*- coding: utf-8 -*-
'''

mnsIZH01.py      FEB 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

IZHIKEVICH MODEL FOR ACTION POTENTIALS AND SPIKE TRAINS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mnsIZH01.pdf

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
def VDOT(t,v,u,Iext):
    Vdot  = c1*v**2 + c2*v + c3 - u + Iext
    return Vdot    
def UDOT(t,v,u,Iext):
    Udot  = a*(b*v - u)
    return Udot   


#%% INPUTS AND MODEL PARAMETERS

a = 0.1
b = 0.2
c = -50
d = 2

v0 = -65
I0 = 10

c1 = 0.04; c2 = 5; c3 = 140

# Time span
N = 9999

tS = 100

t = linspace(0,tS,N)
h = t[2] - t[1]

v = zeros(N); u = zeros(N); Iext = zeros(N)
Iext[t>10] = I0
v[0] = v0
u[0] = b*v[0]

       
#%% Runge-Kutta
for n in range(N-1):
    k1 = VDOT(t[n],v[n],u[n],Iext[n])
    g1 = UDOT(t[n],v[n],u[n],Iext[n])
    
    k2 = VDOT(t[n]+0.5*h,v[n]+k1*h/2,u[n]+g1*h/2,Iext[n])
    g2 = UDOT(t[n]+0.5*h,v[n]+k1*h/2,u[n]+g1*h/2,Iext[n])
    
    k3 = VDOT(t[n]+0.5*h,v[n]+k2*h/2,u[n]+g2*h/2,Iext[n])
    g3 = UDOT(t[n]+0.5*h,v[n]+k2*h/2,u[n]+g2*h/2,Iext[n])
    
    k4 = VDOT(t[n]+1.0*h,v[n]+k3*h,u[n]+g3*h,Iext[n])
    g4 = UDOT(t[n]+1.0*h,v[n]+k3*h,u[n]+g3*h,Iext[n])
    
    v[n+1] = v[n] + h*(k1 + 2*k2 + 2*k3 + k4)/6
    u[n+1] = u[n] + h*(g1 + 2*g2 + 2*g3 + g4)/6
    
    
    if v[n+1] > 30:
    #   V[n] = vSpike
        v[n+1] = c
        u[n+1] = u[n+1] + d
    
#%%   FIG 1: t vs V   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('v  [ mV ]')
ax.set_title('a = %0.2f' %a + '   b = %0.2f' %b +
            '  c = %0.2f' %c + '   d = %0.2f' %d +
            '  I$_0$ = %0.2f' %I0 , fontsize = 11)
ax.grid()
ax.plot(t,v,'b', lw = 2)

fig1.tight_layout()
fig1.savefig('a1.png')

#%%   FIG 2: t vs IEXT   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2)
fig2, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel(r'I$_{ext}$  [ mV.ms$^{-1}$ ]')
ax.grid()
ax.set_ylim([0,11])
ax.plot(t,Iext,'r',lw = 2)

fig2.tight_layout()
fig2.savefig('a2.png')

#%%   FIG 3: t vs u   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('u  [ mV ]')
ax.set_title('a = %0.2f' %a + '   b = %0.2f' %b +
            '  c = %0.2f' %c + '   d = %0.2f' %d +
            '  I$_0$ = %0.0f' %I0 , fontsize = 11)
ax.grid()

ax.plot(t,u,'b', lw = 2)

fig3.tight_layout()
fig3.savefig('a3.png')


#%%   FIG 4: u vs v   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('u  [ Vs ]'); ax.set_ylabel('v  [ mV ]')
ax.set_title('a = %0.2f' %a + '   b = %0.2f' %b +
            '  c = %0.2f' %c + '   d = %0.2f' %d +
            '  I$_0$ = %0.0f' %I0 , fontsize = 11)
ax.grid()

ax.plot(u,v,'b', lw = 2)

fig4.tight_layout()
fig4.savefig('a4.png')

         
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


