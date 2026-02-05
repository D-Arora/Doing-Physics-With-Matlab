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
#plt.close('all')

#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x = state
    m_inf = 1/(1 + np.exp((Vh - x)/k))
    GNa = GNa_max*m_inf
    dx = -( I0 + GL*(x - EL) + GNa*(x - ENa) ) / CM
    return dx  

    
#%% Time span   [s]
N = 9999
tS = 10e-3
t = linspace(0,tS,N)
h = t[2] - t[1]

# Initial membrane voltage  [V]  / ext current [A]
V0 = 100e-3

I0 = -0.02e-3
Iext = ones(N)*I0
#Iext = (I0/tS)*t

# Model parameters
CM = 10e-6      # Membrane capacitance [F]

GL = 19e-3      # leak conductance [S]
GNa_max = 74e-3    # Na+ conductance  [S]

EL = -67e-3     # Leak reversal potential / Nernst potential  [V]
ENa = 60e-3     # Na+ reversal potential
Vh = 19e-3      # m_inf
k = 9e-3        # m_inf


#%% # Solve ODE
sol = odeint(lorenz, V0, t, tfirst=True)
VM = sol[:,0] 

IL = GL*(VM - EL)
m_inf = 1/(1 + np.exp((Vh - VM)/k))
GNa = GNa_max*m_inf
INa = GNa*(VM - ENa)
IC = -(Iext + IL + INa)
Inet = Iext + IL + INa + IC
Vdot = (Iext -GL*(VM - EL) - GNa*(VM - ENa) ) / CM

plt.close('all')
#xP = t*1e3; yP = VM*1e3
#plt.plot(xP,yP,lw = 2)
#XX
    
#%%   FIG 1: t vs VM   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel('V$_M$  [ mV ]')
ax.set_xlabel('t  [ ms ]')
q = I0*1e3;ax.set_title('I$_{ext}$ = %0.2f mA' %q)
ax.grid()

xP = t*1e3; yP = VM*1e3
ax.plot(xP,yP,'b',lw = 2)

#xP = [0,tS*1e3]; yP = [VM[-1]*1e3,VM[-1]*1e3]; ax.plot(xP,yP,'m',lw = 1)
fig1.tight_layout()
fig1.savefig('a1.png')

#xxx

#%%   FIG 2: currents  t vs I 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel('I  [ mA ]')
ax.set_xlabel('t [ ms ]')
ax.grid()
xP = t*1e3; yP = Iext*1e3; ax.plot(xP, yP,'m',lw = 2, label = 'I$_{ext}$')
yP = IL*1e3; ax.plot(xP, yP,'r',lw = 2,label = 'I$_L$')
yP = INa*1e3; ax.plot(xP, yP,'b',lw = 2,label = 'I$_{Na}$')
yP = IC*1e3; ax.plot(xP, yP,'k',lw = 2,label = 'I$_C$')
yP = Inet*1e3; ax.plot(xP, yP,'g',lw = 2,label = 'I$_{net}$')
ax.legend(fontsize = 10, ncols = 5)
fig2.tight_layout()
fig2.savefig('a2.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


