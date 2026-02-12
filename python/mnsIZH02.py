# -*- coding: utf-8 -*-
'''

mnsIZH02.py      FEB 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

IZHIKEVICH MODEL FOR ACTION POTENTIALS AND SPIKE TRAINS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mnsIZH02.pdf

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
    Vdot  = (k*(v - vr)*(v - vt) - u + Iext) /C
    return Vdot    
def UDOT(t,v,u,Iext):
    Udot  = a*(b*(v - vr) - u)
    return Udot   


#%% INPUTS AND MODEL PARAMETERS
I0 = 600
tS = 500

flag = 3
if flag == 1:
   C = 100; vr = -60; vt = -40; k = 0.7 
   a = 0.03; b = -2; c = -50; d = 100; vPeak = 35; 
   
if flag == 2:  #  Intrinsically Bursting (IB) neurons 
   C = 150; vr = -75; vt = -45; k = 1.2 
   a = 0.01; b = 5; c = -56; d = 130; vPeak = 50
   
if flag == 3:  #  Chattering (CH) neurons 
   C = 50; vr = -60; vt = -40; k = 1.5 
   a = 0.03; b = 1; c = -40; d = 150; vPeak = 35  
   

# Initial membrane potential [mV]    
v0 = -40

# Time span
N = 9999
t = linspace(0,tS,N)
h = t[2] - t[1]

v = zeros(N); u = zeros(N); Iext = zeros(N)
v[0] = vr
u[0] = b*v[0]


Iext[t>10] = I0    # step input
#Iext = 0.5*t       # ramp input
       
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
    
    
    if v[n+1] > vPeak:
    #   V[n] = vSpike
        v[n+1] = c
        u[n+1] = u[n+1] + d
    
#%%   FIG 1: t vs V   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('v  [ mV ]')
ax.set_title('C = %0.0f' %C+ '    k = %0.2f' %k +
             '    v$_r$ = %0.0f' %vr + '    v$_t$ = %0.0f' %vt + 
             '\n a = %0.2f' %a + '   b = %0.1f' %b +
            '   c = %0.0f' %c + '    d = %0.0f' %d +
            '   I$_0$ = %0.1f' %I0 , fontsize = 11)
ax.grid()
ax.plot(t,v,'b', lw = 2)

fig1.tight_layout()
fig1.savefig('a1.png')

#%%   FIG 2: t vs IEXT   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2)
fig2, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel(r'I$_{ext}$  [ pA ]')
ax.set_title('I$_0$ = %0.0f  pA' %I0)
ax.grid()
ax.set_ylim([0,1.1*max(Iext)])
ax.plot(t,Iext,'r',lw = 2)

fig2.tight_layout()
fig2.savefig('a2.png')

#%%   FIG 3: t vs u   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('u  [ mV ]')

ax.set_title('C = %0.0f' %C+ '    k = %0.2f' %k +
             '    v$_r$ = %0.0f' %vr + '    v$_t$ = %0.0f' %vt + 
             '\n a = %0.2f' %a + '   b = %0.1f' %b +
            '   c = %0.0f' %c + '    d = %0.0f' %d +
            '   I$_0$ = %0.1f' %I0 , fontsize = 11)
ax.grid()

ax.plot(t,u,'b', lw = 2)

fig3.tight_layout()
fig3.savefig('a3.png')


#%%   FIG 4: u vs v   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('u  [ pA ]'); ax.set_ylabel('v  [ mV ]')

ax.set_title('C = %0.0f' %C+ '    k = %0.2f' %k +
             '    v$_r$ = %0.0f' %vr + '    v$_t$ = %0.0f' %vt + 
             '\n a = %0.2f' %a + '   b = %0.1f' %b +
            '   c = %0.0f' %c + '    d = %0.0f' %d +
            '   I$_0$ = %0.1f' %I0 , fontsize = 11)

ax.grid()

ax.plot(u,v,'b', lw = 1)
ax.plot(u[0],v[0],'bo', ms = 6)

fig4.tight_layout()
fig4.savefig('a4.png')


#%%   Response to a ramp input: ISI 
peaks = find_peaks(v)
Peaks = peaks[0]
L = len(Peaks)
tPeaks = t[Peaks]
tPeaks2 = tPeaks[1:L]
tPeaks1 = tPeaks[0:L-1]
ISI = tPeaks2 - tPeaks1
f = 1000/ISI
tPEAKS = (tPeaks2+tPeaks1)/2
index = (Peaks[0:L-1] + Peaks[1:L])/2
index = index.astype(int)
IPEAKS = Iext[index]

print('I0 = %0.1f' %I0)
q = np.mean(f); print('mean f = %0.2f' %q)
q = np.mean(ISI); print('mean ISI = %0.0f' %q)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig5, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('$I_{ext}$  [ pA ]'); ax.set_ylabel('f  [ Hz ]')

ax.set_title('C = %0.0f' %C+ '    k = %0.2f' %k +
             '    v$_r$ = %0.0f' %vr + '    v$_t$ = %0.0f' %vt + 
             '\n a = %0.2f' %a + '   b = %0.1f' %b +
            '   c = %0.0f' %c + '    d = %0.0f' %d 
            , fontsize = 11)

ax.grid()
ax.set_xlim([0,max(IPEAKS)])
ax.set_ylim([0,1.1*max(f)])

ax.plot(IPEAKS,f,'go', ms = 3)

fig5.tight_layout()
fig5.savefig('a5.png')
         
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


