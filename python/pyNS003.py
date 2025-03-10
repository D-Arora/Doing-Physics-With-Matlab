# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:33:51 2025

@author: Owner
"""

"""
pyNS001.py          Jan 2025

NEUROSCIENCE
ION CHANNEL DYNAMICS: voltage clamp ion conductances and currents
 
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/????.pdf
    

 
"""

#%%
import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array, ones 
import matplotlib.pyplot as plt
from scipy.integrate import odeint, quad, dblquad, simps
import time
from sympy import solve, symbols
from scipy import signal
tStart = time.time()

#%%  INPUTS / SIMULATION TIME
N = 9999    # <<<<<  Grid points

tMax = 50   # <<<<<  Simulation time [ms] 

t = np.linspace(0,tMax,N)
dt = t[2] - t[1]  

#%%   FIXED PARAMETERS 
T = 20            # temperature [20 deg C]
eps = 1e-16

EN = 50          # reversal voltage Na+ [mV]
EK = -77         # reversal voltage K+  [mV]
EL = -76         # reversal voltage leak [mV]
C  = 1.0          # membrane capacitance/area  [uF.cm^-2]

GK = 36          # K+ max conductance [mS.cm^-2]
GN = 120         # Na+ conductance [mS.cm.-2)]
GL = 0.3         # max leakage conductance [mS.cm-2]

#%% SETUP TIME DEPENDENT MATRICES
V    = zeros(N)       # membrane potential (mV)
JN   = zeros(N)       # Na+ current density (uA.cm^-2)
JK   = zeros(N)       # K+  current density (uA.cm^-2)
JL   = zeros(N)       # leakage current density (uA.cm^-2)
Jm   = zeros(N)       # membrane current (uA.cm^-2)
Jext = zeros(N)       # External current stumulus  (uA.cm-2)

gN = zeros(N)       # Na+ conductance
gK = zeros(N)       # K+ conductance
gL = GL*ones(N)        # gL conductance

n    = zeros(N)       # K+ gate parameter
m    = zeros(N)       # Na+ gate parameter
h    = zeros(N)       # Na+ gate parameter

#%% EXTERNAL CURRENT STIMULUS
flag = 4         # <<<<< Select external current 

if flag == 1:    # Pulse   t1 (on) / t2 (off)
   J0 = 20       # Amplitude of pulse 
   t1 = 0.5; t2 = 1.0 
   N1 = round(t1/dt); N2 = round(t2/dt)
   Jext[N1:N2] = J0

if flag == 2:    # Dual pulses
   J0 = 20       # Amplitude of pulse 
   p = 0.2       # pulse width
   t1 = 0.5; t2 = t1+p; 1.0; t3 = 7.0; t4 = t3+p 
   Jext[t>t1] = J0; Jext[t>t2] = 0; Jext[t>t3] = J0; Jext[t>t4] = 0 
   
if flag == 3:    # Square wave: pulse train
   J0 = 100     # Amplitude of pulse 
   p = 0.1        # period [ms]
   Jext = J0*sin(2*pi*t/p)
   Jext[Jext>0] = J0; Jext[Jext<0] = 0   
   
if flag == 4:
    J0 = 30
    t1 = 5
    Jext[t>t1] = J0


     

#%% INITIAL VALUES
phi = 3**((T-6.3)/10)    # temperature dependent variable
dV = 0
An0 = phi * (eps + 0.10 - 0.01 * dV) / (eps + exp(1 - 0.1 * dV) - 1)
Am0 = phi * (eps + 2.5 - 0.1  * dV)  / (eps + exp(2.5 - 0.1 * dV) - 1)
Ah0 = phi * 0.07 * exp(-dV / 20)

Bn0 = phi * 0.125 * exp(-dV / 80)
Bm0 = phi * 4 * exp(-dV/18)
Bh0 = phi * 1 / (exp(3.0 - 0.1 * dV) + 1)

n0 = An0/(An0+Bn0); n[0] = n0
m0 = Am0/(Am0+Bm0); m[0] = m0
h0 = Ah0/(Ah0+Bh0); h[0] = h0

gN0 = GN*m0**3*h0; gN[0] = gN0
gK0 = GK*n0**4   ; gK[0] = gK0
gL0 = gL[0]      

x = symbols('x')   # find resting membrane potential
z = solve(gN0*(x-EN)+gK0*(x-EK)+gL0*(x-EL), x)
Vrest = z[0]
Vrest = float(Vrest)
V[0] = Vrest

JN0 = gN0*(Vrest - EN); JN[0] = JN0
JK0 = gK0*(Vrest - EK); JK[0] = JK0
JL0 = gL0*(Vrest - EL); JL[0] = JL0
Jm0 = JN0 + JK0 + JL0 ; Jm[0] = Jm0

#%% SOLVE ODE
for s in range(N-1):
    
    dv = V[s] - Vrest
    An = phi * (eps + 0.10 - 0.01 * dv) / (eps + exp(1 - 0.1 * dv) - 1)
    Am = phi * (eps + 2.5 - 0.1  * dv)  / (eps + exp(2.5 - 0.1 * dv) - 1)
    Ah = phi * 0.07 * exp(-dv / 20)

    Bn = phi * 0.125 * exp(-dv / 80)
    Bm = phi * 4 * exp(-dv/18)
    Bh = phi * 1 / (exp(3.0 - 0.1 * dv) + 1)

    n[s+1] = n[s] + dt * (An *(1-n[s]) - Bn * n[s]) 
    m[s+1] = m[s] + dt * (Am *(1-m[s]) - Bm * m[s]) 
    h[s+1] = h[s] + dt * (Ah *(1-h[s]) - Bh * h[s]) 

    gK[s+1] = n[s+1]**4 * GK
    gN[s+1] = m[s+1]**3 * h[s+1] * GN

    JK[s+1] = gK[s+1] * (V[s] - EK)
    JN[s+1] = gN[s+1] * (V[s] - EN)
    JL[s+1] = gL[s+1] * (V[s] - EL)
   
    V[s+1] = V[s] + (dt/C) * (-JK[s+1] - JN[s+1] - JL[s+1] + Jext[s+1])
  
Jm = JN + JK + JL
vDot = np.gradient(V,dt)              # dv/dt
vDot2 = np.gradient(vDot,dt)          # d2v/dt2

# Charge delivered by external current Q [ nC.cm-2]
Q = simps(Jext,t)
Vmax = max(V)

#%% CONSOLE OUTPUT
print(' ')
print('Charge delivered by Jext  Q = %0.2f nC' %Q )
print('max membrane potential  Vmax = %0.2f mV' %Vmax )
if flag == 1:
   s = t2-t1; print('Pulse duration t2 - t1 = %0.2f ms' %s ) 
   print('Pulse heigth  J0 = %0.2f uA' %J0 ) 
if flag == 2:
   s = t3-t2; print('Time delay interval between pulses = %0.2f ms' %s ) 
   print('Pulse heigth  J0 = %0.2f uA' %J0 )    
s = time.time() - tStart; print('Execution time = %2.2f' %s) 
  
#%% GRAPHICS
# Fig 1
plt.rcParams["figure.figsize"] = (6,2.4)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t [ms]',fontsize = 12)
ax.set_ylabel('$V_{m}$  [ mV ] ',fontsize = 12)
ax.grid()
xP = t; yP = V
ax.plot(xP,yP,'b',lw = 2)
ax2 = ax.twinx()
xP = t; yP = Jext
ax2.plot(xP,yP,'r',lw = 0.4)
ax2.set_ylabel('$J_{ext}$  $\mu$A.cm$^{-2}$ ',fontsize = 12, color = 'r')
fig1.tight_layout()

# # Fig 2
# plt.rcParams["figure.figsize"] = (6,2.2)
# fig2, ax = plt.subplots(nrows=1, ncols=1)

# ax.set_xlabel('t  [ms]',fontsize = 12)
# ax.set_ylabel('$J_{Na}$  [ $\mu$A.cm$^{-2}$ ] ',fontsize = 12)
# ax.grid()
# xP = t; yP = JN
# ax.plot(xP,yP,'b',lw = 2)
# ax2 = ax.twinx()
# xP = t; yP = Jext
# ax2.plot(xP,yP,'r',lw = 0.4)
# ax2.set_ylabel('$J_{ext}$  $\mu$A.cm$^{-2}$ ',fontsize = 12, color = 'r')
# fig2.tight_layout()
# # Fig 3
# plt.rcParams["figure.figsize"] = (6,2.2)
# fig3, ax = plt.subplots(nrows=1, ncols=1)
# ax.set_xlabel('t  [ms]',fontsize = 12)
# ax.set_ylabel('$J_{K}$  [ $\mu$A.cm$^{-2}$ ] ',fontsize = 12)
# ax.grid()
# xP = t; yP = JK
# ax.plot(xP,yP,'b',lw = 2)
# ax2 = ax.twinx()
# xP = t; yP = Jext
# ax2.plot(xP,yP,'r',lw = 0.4)
# ax2.set_ylabel('$J_{ext}$  $\mu$A.cm$^{-2}$ ',fontsize = 12, color = 'r')
# fig3.tight_layout()
# # Fig 4
# plt.rcParams["figure.figsize"] = (6,2.2)
# fig4, ax = plt.subplots(nrows=1, ncols=1)
# ax.set_xlabel('t  [ms]',fontsize = 12)
# ax.set_ylabel('$J_{L}$  [ $\mu$A.cm$^{-2}$ ] ',fontsize = 12)
# ax.grid()
# xP = t; yP = JL
# ax.plot(xP,yP,'b',lw = 2)
# ax2 = ax.twinx()
# xP = t; yP = Jext
# ax2.plot(xP,yP,'r',lw = 0.4)
# ax2.set_ylabel('$J_{ext}$  $\mu$A.cm$^{-2}$ ',fontsize = 12, color = 'r')
# fig4.tight_layout()
# # Fig 5
# plt.rcParams["figure.figsize"] = (6,2.2)
# fig5, ax = plt.subplots(nrows=1, ncols=1)
# ax.set_xlabel('t  [ms]',fontsize = 12)
# ax.set_ylabel('$J_{m}$  [ $\mu$A.cm$^{-2}$ ] ',fontsize = 12)
# ax.grid()
# xP = t; yP = Jm
# ax.plot(xP,yP,'b',lw = 2)
# ax2 = ax.twinx()
# xP = t; yP = Jext
# ax2.plot(xP,yP,'r',lw = 0.4)
# ax2.set_ylabel('$J_{ext}$  $\mu$A.cm$^{-2}$ ',fontsize = 12, color = 'r')
# fig5.tight_layout()
# # Fig 6
# plt.rcParams["figure.figsize"] = (6,2.2)
# fig6, ax = plt.subplots(nrows=1, ncols=1)
# ax.set_xlabel('t  [ms]',fontsize = 12)
# ax.set_ylabel('$J_{ext}$  [ $\mu$A.cm$^{-2}$ ] ',fontsize = 12)
# ax.grid()
# xP = t; yP = Jext
# ax.plot(xP,yP,'b',lw = 2)
# fig6.tight_layout()
# fig 7
plt.rcParams["figure.figsize"] = (6,3)
fig7, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ms]',fontsize = 12)
ax.set_ylabel('n  m  h ',fontsize = 12)
ax.grid()
xP = t; yP = n
ax.plot(xP,yP,'b',lw = 2,label = 'n')
yP = m
ax.plot(xP,yP,'r',lw = 2,label = 'm')
yP = h
ax.plot(xP,yP,'k',lw = 2,label = 'h')
ax2 = ax.twinx()
xP = t; yP = Jext
ax2.plot(xP,yP,'r',lw = 0.4)
ax2.set_ylabel('$J_{ext}$  $\mu$A.cm$^{-2}$ ',fontsize = 12, color = 'r')
ax.legend()
fig7.tight_layout()
# fig 8
plt.rcParams["figure.figsize"] = (6,3)
fig8, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ms]',fontsize = 12)
ax.set_ylabel('g  [ mS.cm$^{-2}$ ] ',fontsize = 12)
ax.grid()
xP = t; yP = gN
ax.plot(xP,yP,'b',lw = 2,label = 'Na+')
yP = gK
ax.plot(xP,yP,'r',lw = 2,label = 'K+')
ax2 = ax.twinx()
xP = t; yP = Jext
ax2.plot(xP,yP,'r',lw = 0.4)
ax2.set_ylabel('$J_{ext}$  [ $\mu$A.cm$^{-2}$ ] ',fontsize = 12, color = 'r')
ax.legend()
fig8.tight_layout()

# fig 9
plt.rcParams["figure.figsize"] = (6,3)
fig9, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('$V_m$ [mV]',fontsize = 12)
ax.set_ylabel('dV/dt [ mV.ms$^{-1}$ ] ',fontsize = 12)
ax.grid()
xP = V; yP = vDot
ax.plot(xP,yP,'b',lw = 2)
xP = V[0]; yP = vDot[0]
ax.plot(xP,yP,'og',ms = 8)
xP = V[3000]; yP = vDot[3000]
ax.plot(xP,yP,'om',ms = 8)
fig9.tight_layout()

# fig 10 ion currents
plt.rcParams["figure.figsize"] = (6,4)
fig10, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ms]',fontsize = 12)
ax.set_ylabel('$J_{m}$  [ $\mu$A.cm$^{-2}$ ] ',fontsize = 12)
ax.grid()
xP = t; yP = Jm
ax.plot(xP,yP,'k',lw = 2, label = 'J$_m$')
xP = t; yP = JN
ax.plot(xP,yP,'b',lw = 2, label = 'J$_{NA}$')
xP = t; yP = JK
ax.plot(xP,yP,'r',lw = 2, label = 'J$_K$')
xP = t; yP = JL
ax.plot(xP,yP,'m',lw = 2, label = 'J$_{L}$')

ax2 = ax.twinx()
xP = t; yP = Jext
ax2.plot(xP,yP,'r',lw = 0.4)
ax2.set_ylabel('$J_{ext}$  $\mu$A.cm$^{-2}$ ',fontsize = 12, color = 'r')
ax.legend()
fig10.tight_layout()


#%% SAVE PLOTS
fig1.savefig('a1.png')
# fig2.savefig('a2.png')
# fig3.savefig('a3.png')
#fig4.savefig('a4.png')
# fig5.savefig('a5.png')
# fig6.savefig('a6.png')
fig7.savefig('a7.png')
fig8.savefig('a8.png')
fig9.savefig('a9.png')
fig10.savefig('a10.png')





