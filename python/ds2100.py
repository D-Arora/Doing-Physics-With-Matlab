# -*- coding: utf-8 -*-
'''
ds2100.py      Oct 2025
DYNAMICAL SYSTEMS
LIMIT CYCLES: VAN DER POL OSCILLATOR

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS.ds2100.pdf

Displacement  x
Velocity      v = xDot = y

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time

tStart = time.time()
plt.close('all')

#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x, y = state
    dx = y
    dy = -x + mu*(1 - x**2)*y
    return [dx, dy]  

#%%  INPUTS >>>
# Damping coefficient
mu = 0
# Initial conditions xI yI (vI)
xI,yI = 2,0
# time span
tS = 60; nT = 9999
# Plot color
col = [0,0,1]

#%% Solution ODE for x and y 
t = linspace(0,tS,nT)
u0 = [xI,yI]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]       

peaks, _ = find_peaks(xS, height = 1.5)
T = t[peaks[-1]] - t[peaks[-2]]
print('period  T = %0.2f' %T)

# Steady sate amplitudes
xSS = max(xS[t>0.50*tS])
ySS = max(yS[t>0.50*tS])

#%% PHASE PORTRAIT streamplot
N = 11
x = linspace(-4,4,N)
xx,yy = np.meshgrid(x,x)
xxDot = yy
yyDot = -xx + mu*(1 - xx**2)*yy

#%%   FIG 1: t vs x   t vs y = v
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,3.5)
fig1, ax = plt.subplots(nrows=2, ncols=1)

C = 0   
ax[C].set_xlabel('t'); ax[C].set_ylabel('x')
ax[C].grid()
ax[C].set_title('$\mu$ = %0.2f' %mu + '   x$_{SS}$ = %0.2f' %xSS)
ax[C].plot(t,xS,lw = 2,color = col)

C = 1   
ax[C].set_xlabel('t'); ax[C].set_ylabel('v')
ax[C].grid()
ax[C].set_title('$\mu$ = %0.2f' %mu + '  v$_{SS}$ = %0.2f' % ySS)
ax[C].plot(t,yS,lw = 2,color = col)

fig1.tight_layout()
fig1.savefig('a1.png')

#%% FIGURE 2: Phase Portrait  quiver
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,3)
fig2, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x')
axes.set_ylabel('v', rotation = 0)
axes.set_title('$\mu$ = %0.2f' %mu )
axes.grid()

axes.quiver(xx,yy,xxDot,yyDot,linewidth = 1)
axes.plot(xS,yS,lw = 2,color = col)
axes.plot(0,0,'mo',ms = 6)
axes.plot(xS[0],yS[0],'go',ms = 6)

fig2.tight_layout()
fig2.savefig('a2.png')

#%% FIGURE 3: Phase Portrait  streamplot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3,2) #(4,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x')
axes.set_ylabel('v', rotation = 0)
axes.set_title('$\mu$ = %0.2f' %mu + '   T = %0.1f' %T )
axes.grid()

#axes.streamplot(xx,yy,xxDot,yyDot,linewidth = 1, density = 1.2)
axes.plot(0,0,'mo',ms = 6)

fig3.tight_layout()
fig3.savefig('a3.png')

#xxx
#%% SOLVE ODE
# Time span 
t1 = 0; t2 = 60; nT = 9999
t = linspace(t1,t2,nT)
N = 1    # number of initial condition points
x0 = zeros(N); y0 = zeros(N)

x0[0], y0[0] = xI,yI
#x0[1], y0[1] = -2,10
#x0[2], y0[2] = 2,-10
# x0[3], y0[3] = 1.5,0
# x0[4], y0[4] = 2,0

for c in range(N):
  u0 = [x0[c],y0[c]]
  sol = odeint(lorenz, u0, t, tfirst=True)
  xS = sol[:,0]     
  yS = sol[:,1]       
  
  plt.plot(xS,yS,lw = 2)
  plt.plot(xS[0],yS[0],'go',ms = 6) 
    
fig3.savefig('a3A.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


