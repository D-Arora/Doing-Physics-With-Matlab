# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 21:22:37 2024

@author: Owner
"""

# https://math.libretexts.org/Bookshelves/Scientific_Computing_Simulations_and_Modeling/Scientific_Computing_(Chasnov)/II%3A_Dynamical_Systems_and_Chaos/12%3A_Concepts_and_Tools

# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace, zeros

plt.close('all')
tStart = time.time()

#%%  SOLVE ODE
# Solve ODE for x
def lorenz(t, state):    
    xS = state
    dx = r*xS + xS**3 - xS**5
    return dx  

#%%  Figure 1: x vs xDot  
# Input r value / plot limits>>>
r = -0.3

y1 = -0.3; y2 = 0.3
# Zeros of xDot 
coeff = [-1,0,1,0,r]     
Z = np.roots(coeff)

num = 99
x1 = -1.2; x2 = 1.2; x = linspace(x1,x2,num)
xDot = r*x + x**3 - x**5

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)

fig1, ax = plt.subplots(1)
ax.set_xlabel('x',color = 'black')
ax.set_ylabel('$x_{dot}$    f(x)',color= 'black') 
ax.set_title('r = %2.2f' %r ) 
 
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_xlim([x1, x2])
ax.set_ylim([y1,y2])
#ax.set_xticks(np.arange(-4,4,1))

xP = x; yP = xDot; 
ax.plot(xP,yP,'k',lw = 2, label = 'r = -10')

xP = Z; yP = zeros(len(Z))
if r > -0.25: ax.plot(xP,yP,'ko',ms = 6)
ax.plot(0,0,'ko',ms = 6)
fig1.tight_layout()

fig1.savefig('a1.png')

#%% Figure 2: time evolution t vs x
# INPUTS r and x(0)  >>>>>
r = 0
x0 = 0.05

y1 = -20; y2 = 0; dy = 5

N = 9999; t1 = 0; t2 = 300; t = linspace(t1,t2,N)

sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,2.8)
fig2, ax = plt.subplots(1)
ax.set_xlabel('t',color = 'black')
ax.set_ylabel('x',color= 'black')
ax.set_title(r'r = %2.0f' %r  + '   x(0) = %2.5f' %x0, fontsize = 12)
ax.grid()

#ax.set_ylim([-1, 1])
#ax.set_yticks(np.arange(y1,y2,dy))

xP = t; yP = xS; 
ax.plot(xP,yP,'b',lw = 2)    

fig2.tight_layout()
fig2.savefig('a2.png')

#%% Figure 3: Bifurcation diagram

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(1)

N = 999
r1 = -0.5; r2 = 0
xP = linspace(r1,r2,N); yP = zeros(N)
ax.plot(xP,yP,'b',lw = 2) 

r1 = 0; r2 = 1
xP = linspace(r1,r2,N); yP = zeros(N)
ax.plot(xP,yP,'r',lw = 2) 

r1 = linspace(-1/4,1,999)
z1 = 1 + (1 + 4*r1)**0.5
xe1 = (0.5*z1)**0.5
xP = r1; yP = xe1
ax.plot(xP,yP,'b',lw = 2) 
ax.plot(xP,-yP,'b',lw = 2) 

r2 = linspace(-1/4,0,999)
z2 = 1 - (1 + 4*r2)**0.5
xe3 = (0.5*z2)**0.5
xe4 = -xe3
xP = r2; yP = xe3
ax.plot(xP,yP,'r',lw = 2) 
yP = xe4
ax.plot(xP,yP,'r',lw = 2) 

ax.set_xlabel('r',color = 'black',fontsize = 12)
ax.set_ylabel('$x_e$',color= 'black',fontsize = 14)
ax.grid()
ax.set_xticks(np.arange(-0.5,1.1,0.25))

fig3.tight_layout()      
fig3.savefig('a3.png')
