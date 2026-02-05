# -*- coding: utf-8 -*-
'''
mns007.py      Jan 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE
   LEAK / FAST Na+ Channels
GEOMETRICAL ANAYLSIS OF [1D] DYNAMICAL SYSTEMS

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
# Iext == I0
I0 = -0.02e-3

# Model parameters
CM = 10e-6      # Membrane capacitance [F]
GL = 19e-3      # leak conductance [S]
GNa_max = 74e-3    # Na+ conductance  [S]
EL  = -67e-3     # Leak reversal potential / Nernst potential  [V]
ENa = 60e-3      # Na+ reversal potential
Vh = 19e-3       # m_inf
k  = 9e-3        # m_inf

#%%
Vss = (-I0 + GNa_max*ENa + GL*EL) / (GNa_max + GL)

#%% Membrane potential VM == V
V1 = -100e-3; V2 = 60e-3
V = linspace(V1,V2,N)
m_inf = 1/(1 + np.exp((Vh - V)/k))
GNa = GNa_max * m_inf
# vDot  [V/s]
Vdot = -(I0 + GL*(V - EL) + GNa*(V - ENa) ) / CM

# Find V for Vdot = 0
z = zeros(3); p = 0; 
for c in range(N-2):
    q = Vdot[c]*Vdot[c+1]
    if q <= 0:
       z[p] = c
       p = int(p+1)     
z = z.astype(int)
Vzero = V[z]*1e3
print(np.round(Vzero,2))

Iext = I0*ones(N)  
IL = GL*(V - EL)
INa = GNa*(V - ENa)
IC = -(I0 + IL + INa)
Inet = Iext + IL + INa + IC 


#%%   FIG 1: Phae plot  vM vs Vdot   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel('V$_{dot}$  [ mV/ms ]')
ax.set_xlabel('V$_M$  [ mV ]')
q = I0*1e3;ax.set_title('I$_{ext}$ = %0.2f mA' %q)
ax.grid()
ax.plot(V*1e3,Vdot,'b',lw = 2)
xP = array([V1,V2])*1e3; yP = [0,0]
ax.plot(xP,yP,'r',lw = 1)

for c in range(3):
    if Vzero[c] > -100:
       ax.plot(Vzero[c],0,'go',ms = 4)

fig1.tight_layout()
fig1.savefig('a1.png')


#%%   FIG 2: Sodium conductance  
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel('G$_{Na}$  [ mS ]')
ax.set_xlabel('V$_M$  [ mV ]')
ax.grid()
xP = V*1e3; yP = GNa*1e3; ax.plot(xP,yP,'b',lw = 2)
xP = array([Vh,Vh])*1e3; yP = array([0,max(GNa)/2])*1e3
ax.plot(xP,yP,'r',lw = 1)
xP = array([V1,Vh])*1e3; yP = array([max(GNa)/2,max(GNa)/2])*1e3
ax.plot(xP,yP,'r',lw = 1)
fig2.tight_layout()
fig2.savefig('a2.png')

#%%   FIG 3: I-V characteristic I vs V 

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel(r'I  [ mA ]')
ax.set_xlabel('V$_M$ [ mV ]')
q = I0*1e3;ax.set_title('I$_{ext}$ = %0.2f mA' %q)
ax.grid()
#ax.set_ylim([-2,5])

xP = V*1e3; yP = IL*1e3; ax.plot(xP, yP,'r',lw = 2, label = 'I$_L$')
yP = INa*1e3; ax.plot(xP, yP,'b',lw = 2, label = 'I$_{Na}$')
yP = IC*1e3; ax.plot(xP, yP,'m',lw = 2, label = 'I$_C$')
yP = Iext*1e3*ones(N); ax.plot(xP, yP,'k',lw = 2, label = 'I$_{ext}$')
yP = Inet*1e3; ax.plot(xP, yP,'g',lw = 1, label = 'I$_{net}$')

for c in range(3):
    if Vzero[c] > -100:
       ax.plot(Vzero[c],0,'ko',ms = 6)

ax.legend(frameon=False,fontsize = 10, ncols = 5)
fig3.tight_layout()
fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


