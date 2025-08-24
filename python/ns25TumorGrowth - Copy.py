# -*- coding: utf-8 -*-
"""
nsTumorGrowth.py

NONLINEAR [1D] DYNAMICAL SYSTEMS
MODELLING TUMOR GROWTH AND TREATMENT


# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/ns25TumorGrowth.pdf

"""

# LIBRARIES  ================================================================
import numpy as np
from numpy import pi, linspace
from scipy.integrate import odeint
from scipy import integrate

import matplotlib.pyplot as plt
import time

tStart = time.time()
plt.close('all')


#%%  SOLVE ODE 
def lorenz(t, state):
    V = 0
    if t > 20: V = V0 
    if t > 30: V = 0
    # if t > 20: V = V0
    # if t > 25: V = 0
    # if t > 30: V = V0
    # if t > 35: V = 0
    # if t > 40: V = V0
    # if t > 45: V = 0
    # if t > 50: V = V0
    # if t > 55: V = 0
    
    N, T, I, C = state
    #if T < 1e-6: T = 0
    dN = r2*N*(1-N/K2) - c4*T*N - a3*C*N
    dT = r1*T*(1-T/K1) - c2*I*T - c3*T*N - a2*C*T
    dI = s + k1*I*T/(k2+T) - c1*I*T - d1*I - a1*C*I
    dC = k3*V - d2*C
    return dN, dT, dI, dC  

#%% Model default parameters (DO NOT CHANGE)
r2,r1 = 1, 1.5
K2, K1 = 1,1
a1, a2, a3 = 0.2,0.3,0.1
c1,c2,c3,c4 = 1, 0.5, 1, 1
d1, d2 = 0.2,1.0
s = 0.33
k1, k2, k3 = 0.01, 0.3, 1 
V0 = 0
# Initial conditions  N(0)  T(0) I(0)  
u0 = np.zeros(5); u0 = [1,0.1,0.1,0]
# Time span
num = 9999; tMax = 200 

#%% Change value of model parameters: comment / uncomment
r1 = 1.5
c2 = 0.47
#u0[0] = 0.4
#u0[1] = 0.5
#K2 = 1.0; K1 = 0.2#
V0 = 0
tMax = 200

#%% Solve ode
t = linspace(0,tMax,num)
sol = odeint(lorenz, u0, t, tfirst=True)
Ns = sol[:,0] 
Ts = sol[:,1] 
Is = sol[:,2] 
Cs = sol[:,3]

# Total drug delivered
A = integrate.simpson(Cs,t)

#%% CONSOLE OUTPUT
# Final values  /  volumes  /  diameters
Nf = Ns[-1]; Tf = Ts[-1]; If = Is[-1]; Cf = Cs[-1]
vol = 200; p = 1/3; dia = 2*(3*200/(4*pi))**p
volN = vol*Ns[-1]; volT = vol*Ts[-1]; volI = vol*Is[-1]
diaN = dia*Ns[-1]**p; diaT = dia*Ts[-1]**p; diaI = dia*Is[-1]**p

print('N(0) = %0.2f' %u0[0] + '  T(0) = %0.2f' %u0[1]
      + '  I(0) = %0.2f' % u0[2] + '  C(0) = %0.2f' %u0[3])
print('Nf = %0.2f' % Nf + '  Tf = %0.2f' % Tf 
      + '  If = %0.2f' % If + '  Cf = %0.2f' % Cf) 
print('volN = %0.2f' % volN + '  volT = %0.2f' % volT + '  volI = %0.2f' % volI)
print('diaN = %0.2f' % diaN + '  diaT = %0.2f' % diaT + '  diaI = %0.2f' % diaI)  
print('Drug delivered = %0.0f' %A)    
       
#%%
#FIGURE 1: 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('cells')
ax.grid()
L = 1.1*np.amax([Ns,Ts,Is])
ax.set_ylim([0,L])
#ax.set_yticks([0,0.2,0.4,0.6,0.8,1])
ax.plot(t,Ns,'b',lw = 2, label = 'N')
ax.plot(t,Ts,'r',lw = 2, label = 'T')
ax.plot(t,Is,'m',lw = 2, label = 'I')
ax.legend()
fig1.tight_layout()

#%% FIGURE 2
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.2)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('treatment')
ax.grid()

ax.plot(t,Cs,'r',lw = 2)
fig2.tight_layout()

#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)



