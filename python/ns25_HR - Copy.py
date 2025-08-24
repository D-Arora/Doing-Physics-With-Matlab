# -*- coding: utf-8 -*-
'''
AA\ns\pyCode\ns25_HR.py      July 2025
# COMPLEX SYSTEMS: NEUROSCIENCE
# HINDMARSH-ROSE MODEL BURSTING NEURONS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/ns25_HR.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, real, imag, sqrt 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import random
from scipy.optimize import fsolve
from sympy import symbols, Eq, solve

tStart = time.time()

plt.close('all')


#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x, y, z = state
    I = 0
    if t > T1: I = I0
    if t > T2: I = 0
    dx = y + 3*x**2 - x**3 - z + I
    dy = 1 - 5*x**2 - y
    dz = r*(s*(x-x1) - z)
    return [dx, dy, dz]  


#%% INPUTS
r = 0.005  #0.001
s = 4
x1 = 0.1; y1 = -1.0; z1 = 0.2      # Initial conditions
I0 = 1.5                      # external current stimulus; max
T1 = 5; T2 = 2000             # external current stimulus; off / on times
t1 = 0; t2 = 2000; nT = 9999     # time span for solving ODEs
t = linspace(t1,t2,nT)         # time span for solving ODEs
I = zeros(nT); I[t>T1] = I0; I[t>T2] = 0  # ext. current


#%%   CRITICAL POINTS  xC and yC     z = 0, I = 0
x, y, z = symbols('x y z')
SS = solve([y + 3*x**2 - x**3,1 - 5*x**2 - y], x,y)
SS0 = SS[0]
SS1 = SS[1]
SS2 = SS[2]
xC0 = float(SS0[0]); yC0 = float(SS0[1])
xC1 = float(SS1[0]); yC1 = float(SS1[1])
xC2 = float(SS2[0]); yC2 = float(SS2[1])
xC = array([xC0,xC1,xC2])
yC = array([yC0,yC1,yC2])


#%% Solve ODEs 
xc = -0.5*(1+sqrt(5)) 
u0 = [xc,y1,z1]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1] 
zS = sol[:,2]  

#%%  PHASE PORTRAIT  x vs y   /   nullclines      I0 = constant  z = 0
X = linspace(-2,2,20); Y = linspace(-20,3,20)
xx,yy = np.meshgrid(X,Y)
xxDot = yy + 3*xx**2 - xx**3 + I0 
yyDot = 1 - 5**xx**2 - yy

XXDot = xxDot/(sqrt(xxDot**2 + yyDot**2))
YYDot = yyDot/(sqrt(xxDot**2 + yyDot**2))

# Nullclines   I0 = constant    z = 0
xN = linspace(-2,2,999)
yNx = -3*xN**2 + xN**3 - I0
yNy = 1 - 5*xN**2


#%% GRAPHICS
# FIGURE 1: Membrane potential x vs time t
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('time',color= 'black',fontsize = 12)
ax.set_ylabel('membrane potential x',color = 'black',fontsize = 12)
ax.grid()
ax.plot(t,xS,'b',lw = 2)
ax1 = ax.twinx()
ax1.set_ylabel('I$_{ext}$', color = 'r',fontsize = 12)
ax1.tick_params(axis='y', colors='r')
ax1.plot(t,I,'r',lw = 2)
fig1.tight_layout()

# FIGURE 2: recovery var. y & adaption var. z vs time
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('time',color= 'black',fontsize = 12)
ax.set_ylabel('recovery var. y',color = 'black',fontsize = 12)
ax.grid()
ax.plot(t,yS,'b',lw = 2)
ax1 = ax.twinx()
ax1.set_ylabel('adaption var. z', color = 'r',fontsize = 12)
ax1.tick_params(axis='y', colors='r')
ax1.plot(t,zS,'r',lw = 2)
fig2.tight_layout()


#%% Figure 3: Phase Portrait: streamplot  


plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('membrane pot. x',color= 'black',fontsize = 12)
ax.set_ylabel('recovery var. y',color = 'black',fontsize = 12)
ax.set_xlim([-2,2])
ax.set_ylim([-20,3])

ax.streamplot(xx,yy,xxDot,yyDot)
ax.plot(xN,yNx,'b',lw = 2)
ax.plot(xN,yNy,'r',lw = 2)

# Critical points
for c in range(3):
    ax.plot(xC[c],yC[c],'ko',ms = 6)

# Initial conditions and trajectories
#x1, y1 = 0,-10    
u0 = [x1,y1,z1]
t1 = 0; t2 = 520
t = linspace(t1,t2,nT)
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1] 
zS = sol[:,2]      

ax.plot(xS,yS,'g',lw = 2)
ax.plot(xS[0],yS[0],'go',ms = 7)

fig3.tight_layout()   


#%% Figure 4: Phase Portrait: quiver plot 
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('membrane pot.  x',color= 'black',fontsize = 12)
ax.set_ylabel('recovery var.  y',color = 'black',fontsize = 12)
ax.set_xlim([-2,2])
ax.set_ylim([-20,3])

ax.quiver(xx,yy,XXDot,YYDot)
ax.plot(xN,yNx,'b',lw = 2)
ax.plot(xN,yNy,'r',lw = 2)

for c in range(3):
    ax.plot(xC[c],yC[c],'ko',ms = 6)

fig4.tight_layout()   

#%% Console output
print(' ')
print('I0 = %0.3f' % I0)
print('Critical points')
for c in range(3):
    print('  ( %0.3f' %xC[c] + ',  %0.3f )' %yC[c] )   


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


#%%
def fn(a,b):
    da = b + 3*a**2 - a**3 
    db = 1 - 5*a**2 - b
    return da,db

da,db = fn(0,-10)    
print(da,db)

#%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')



'''
#