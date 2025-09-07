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
from numpy import pi, sin, cos, linspace, zeros
from scipy.integrate import odeint
from scipy.integrate import odeint
from scipy import integrate
import matplotlib.pyplot as plt
import time

plt.close('all')


#%%  SOLVE ODE 
def lorenz(t, state):
    V = V0
    if t > 20: V = V0  
    # if t > 15: V = 0
    if t > 30: V = 0
    # if t > 25: V = 0
    # if t > 30: V = V0
    # if t > 35: V = 0
    # if t > 40: V = V0
    # if t > 45: V = 0
    N, T, I, C = state
    dN = r2*N*(1-N/K2) - c4*T*N - a3*C*N
    dT = r1*T*(1-T/K1) - c2*I*T - c3*T*N - a2*C*T
    dI = s + k1*I*T/(k2+T) - c1*I*T - d1*I - a1*C*I
    dC = k3*V - d2*C
    return dN, dT, dI, dC 
 

#%% Model parameters
r2,r1 = 1, 1.5
K2, K1 = 1,1
a1, a2, a3 = 0.2,0.3,0.1
c1,c2,c3,c4 = 1, 0.5, 1, 1
d1, d2 = 0.2,1.0
s = 0.33
k1, k2, k3 = 0.01, 0.3, 1 

V0 = 4

#c2 = 0.47

# Initial conditions  N(0)  T(0) I(0)  
u0 = np.zeros(5); u0 = [1,0.1,0.1,0]

# Time span
num = 9999; tMax = 200; t = linspace(0,tMax,num)


#%% Figure 1: Solve ode and plot steady-state vs r 
num = 299
X = linspace(0.5,5,num)
Nss = zeros(num); Tss= zeros(num); Iss = zeros(num); Css = zeros(num)
for c in range(num):
    r1 = X[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    Ns = sol[:,0] 
    Ts = sol[:,1] 
    Is = sol[:,2] 
    Cs = sol[:,3]
    
    Nss[c] = Ns[-1]
    Tss[c] = Ts[-1]
    Iss[c] = Is[-1]

# Drug delivered
drug = integrate.simpson(Cs,t)    
 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r$_1$')
ax.set_ylabel('steady-state')
ax.grid()

#ax.set_ylim([0,1.1*max(T)])
ax.plot(X,Nss,'b',lw = 2, label = 'N$_{ss}$')
ax.plot(X,Tss,'r',lw = 3, label = 'T$_{ss}$')



ax.plot(X,Iss,'m',lw = 2, label = 'I$_{ss}$')
ax.legend(ncol=3,loc = 'upper right',fontsize = 10)
fig1.tight_layout()

fig1.savefig('a1.png')

#%%   Figure 2: solve ODE and plot steady-state vs c2 
r1 = 1.5

# Time span
num = 9999; tMax = 200; t = linspace(0,tMax,num)

num = 299
X = linspace(0,1,num)
Nss = zeros(num); Tss= zeros(num); Iss = zeros(num); Css = zeros(num)
for c in range(num):
    c2 = X[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    Ns = sol[:,0] 
    Ts = sol[:,1] 
    Is = sol[:,2] 
    Cs = sol[:,3]
    
    Nss[c] = Ns[-1]
    Tss[c] = Ts[-1]
    Iss[c] = Is[-1]

# Drug delivered
drug = integrate.simpson(Cs,t)    


plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('c$_2$')
ax.set_ylabel('steady-state')
ax.grid()

#ax.set_ylim([0,1.1*max(T)])
ax.plot(X,Nss,'b',lw = 2, label = 'N$_{ss}$')
ax.plot(X,Tss,'r',lw = 3, label = 'T$_{ss}$')
ax.plot(X,Iss,'m',lw = 2, label = 'I$_{ss}$')
#ax.plot(X,Css,'k',lw = 1,label = 'C$_{ss}$')
ax.legend(ncol=3,loc = 'upper left',fontsize = 10)
fig2.tight_layout()

fig2.savefig('a2.png')

#%% FIGURE 3:   time evolution of C
plt.rcParams['font.size'] = 12  
plt.rcParams["figure.figsize"] = (4,2)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('C')
ax.set_title('V$_0$ = %0.1f' %V0 + '  drug = %0.1f' % drug, fontsize = 10)
#ax.set_ylim([0,1.2*max(Cs)])
ax.grid()

#ax.set_ylim([0,1.1*max(T)])
ax.plot(t,Cs,'b',lw = 2)

fig3.tight_layout()

fig3.savefig('a3.png')

