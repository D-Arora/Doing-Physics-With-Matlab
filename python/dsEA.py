# -*- coding: utf-8 -*-
"""
dsE002A.py      DEC 2025
DYNAMICAL SYSTEMS
 MATHEMATICAL EPIDEMILOGY
 mSIR MODELLING

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/dsE002.pdf
"""

import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt, real, imag 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, simpson
import time
import sys

tStart = time.time()
plt.close('all')


#%% INPUT PARAMETERS
# Time grid [day]
nT = 9999
tS = 100

# Total population, initial number of infected and recovered individuals
No = 100
Io = 0.01 
Eo, Ro, Do = 0, 0, 0

# Probability of  S --> E on contact
p = 0.95  # [  ]
# Contact rate S <---> I   [1/day]
c = 2
# Incubation period [day]
Te = 5
# Infectious period [day]
Tb = 10
# Rate R --> S  [day]
s = 0
# Death rate [1/day]
d = 0

#%% Computations
# Rate S--E
a = c*p
# Rate E --> I [1/day]
e = 1/Te
# Rate I --> R [1/day]
b = 1/Tb
# Populations
# Everyone else is susceptible
So = No - Io - Ro - Eo - Do
# Time span
t = np.linspace(0, tS, nT)
dt = t[2] - t[1]

#%% Solve ODEs Euler method
S = zeros(nT); E = zeros(nT); I = zeros(nT); R = zeros(nT)
D = zeros(nT); N = zeros(nT)
S[0] = So; E[0] = Eo; I[0] = Io; R[0] = Ro; D[0] = Do; N[0] = No

for n in range(nT-1):
    #if n == 1000: c = 1; a = c*p 
    N[n] = S[n] + E[n] + I[n] + R[n] - D[n]
    S[n+1] = S[n] + (-a*S[n]*I[n]/N[n] + s*R[n])*dt
    E[n+1] = E[n] + (a*S[n]*I[n]/N[n] - e*E[n])*dt
    I[n+1] = I[n] + (e*E[n] - b*I[n] - d*I[n])*dt
    R[n+1] = R[n] + (b*I[n] - s*R[n])*dt
    D[n+1] = D[n] + (d*I[n])*dt
    
    #if n > 5000: p = 0.95; a = c*p
        
    # if n == 2000:
    #     S[n+1] = S[n+1] + 50
    #     I[n+1] = I[n+1] + 1
   
#%% SETUP
# total live population 
N = S+E+I+R-D
# Basic reproduction number, effective reproduction number
R0 = a/b
Re = R0*(S/No)
# Susceptiblity threshold
ST = S[I==max(I)] 
tST = t[I==max(I)]
# Infectin load
L = simpson(I,t)



#%% Console output
print('   ')
print('p = %0.2f ' % p + '   c = %0.2f ' % c + '  a = %0.2f' % a)
print('Te = %0.2f' % Te + '   e = %0.2f' % e)
print('Tb = %0.2f' % Tb + '   b = %0.2f' % b)
print('d = %0.2f' % d + '   s = %0.2f' % s)
print('R0 = %0.2f' %R0)
print('Initial populations ')
print('  So = %0.2f' % S[0] + '   Eo = %0.2f' %E [0] +
    '   Io = %0.2f' %I[0] + '   Ro = %0.2f' % R[0] +'   Do = %0.2f' % D[0])
print('Final populations ')
print('  Sf = %0.1f' % S[-1] + '   Ef = %0.1f' %E [-1] +
    '   If = %0.1f' %I[-1] + '   Rf = %0.1f' % R[-1] +'   Df = %0.1f' % D[-1])
print('Max Populations')
print('  Smax = %0.1f' % max(S) + '   Emax = %0.1f' %max(E)
      + '   Imax = %0.1f' %max(I) + '   Rmax = %0.1f' %max(R)
      + '   Dmax = %0.1f' %max(D))
print('Peak: ST = %0.1f ' % ST[0] + '   tST = %0.1f days' % tST[0])
print('Infection load: L = %0.0f ' % L) 

#%% FIG 1:  Population dynamics
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,2.6)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [days]'); ax.set_ylabel('% individuals')
ax.set_title('R$_0$ = %0.2f' %R0 )
ax.grid()
xP = t
yP = (S/No)*100; ax.plot(xP,yP,'b',lw = 2,label = 'S')
yP = (I/No)*100; ax.plot(xP,yP,'m',lw = 2,label = 'I')
yP = (R/No)*100; ax.plot(xP,yP,'g',lw = 2,label = 'R')
#yP = (D/No)*100; ax.plot(xP,yP,'r',lw = 2,label = 'D')
yP = (N/No)*100; ax.plot(xP,yP,'k',lw = 1,label = 'N')
yP = (E/No)*100; ax.plot(xP,yP,'y',lw = 2,label = 'E')
#xP = tST; yP = 100*ST/No; ax.plot(xP,yP,'bo',ms = 4)

ax.legend(fontsize = 10,bbox_to_anchor=(1.2, 1), loc='upper right')
fig1.tight_layout()
fig1.savefig('a1.png')


#%% FIG 2:  INFECTIONS, REPRODUCTION NUMBER 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.6)
fig2, ax = plt.subplots(nrows=1, ncols=1)
# ax2 = ax.twinx()
ax.set_xlabel('t  [days]')
# ax2.set_ylabel('R$_e$', fontsize = 14)
# ax.set_xlim([0,tS])
# ax2.grid()
# ax2.yaxis.label.set_color('b') 
# ax2.tick_params(axis='y', colors='b')  
# ax2.plot(t,Re,'b',lw = 2)
# ax2.plot([0, tS],[1,1],'g',lw = 1.6)

#ax.yaxis.label.set_color('m') 
#ax.tick_params(axis='y', colors='m')  
ax.set_ylabel(' % I ',fontsize = 14)

ax.plot(t,100*I/No,'m',lw = 2)

ax.grid()
fig2.tight_layout()
fig2.savefig('a2.png')

#%%  FIG 3:  phase portrait S vs I
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('S'); ax.set_ylabel('I')
ax.grid()
    
xP = 100*S/No; yP = 100*I/No    
ax.plot(xP,yP,'b',lw = 2,label = 'R')
ax.plot(xP[0],yP[0],'go',ms = 6)

fig3.tight_layout()
# fig3.savefig('a3.png')



#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)