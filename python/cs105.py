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
from numpy import pi, sin, cos, linspace




tStart = time.time()


#%%   FUNCTIONS
def funct(r,x):
    xDot = r*x + x**3 #- x**5
   
    return xDot

#%%  SOLVE ODE
# Solve ODE for x
def lorenz(t, state):    
    xS = state
    dx = r*xS + xS**3 - xS**5
    return dx  


r = -0.25

#%%  CELL 1   
col = [0,0,1]
x1 = -2
x2 = 2
x = linspace(x1,x2,1999)

xDot = funct(r,x)

# if r == 0:
#    xe = 0
# if r < 0:
#    xe = r

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)

fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.20,\
                    right = 0.92, hspace = 0.20,wspace=0.2)

ax.set_xlabel('x',color = 'black')
ax.set_ylabel('$x_{dot}$    f(x)',color= 'black') 
ax.set_title('r = %2.2f' %r ) 
 
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_xlim([x1, x2])
ax.set_ylim([-10, 10])
#ax.set_xticks(np.arange(-4,4,1))

xP = x; yP = xDot; 
ax.plot(xP,yP,color = col,lw = 2, label = 'r = -10')

# xP = -3; yP = 0; col1 = [1,0,0]
# ax.plot(xP,yP,'o',color = col1, ms = 6)
# xP = 0; yP = 0; col1 = [0,0,1]
# ax.plot(xP,yP,'o',color = col1, ms = 6)
# xP = 3; yP = 0; col1 = [1,0,0]
# ax.plot(xP,yP,'o',color = col1, ms = 6)

fig.savefig('a1.png')


#%% CELL 2
r = 0
x0 = 0.05
y1 = -20
y2 = 0
dy = 5

N = 9999
t1 = 0
t2 = 300
t = linspace(t1,t2,N)

sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,2.8)

fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.90, bottom = 0.23, left = 0.22,\
                    right = 0.95, hspace = 0.20,wspace=0.2)
  
ax.set_xlabel('t',color = 'black')
ax.set_ylabel('x',color= 'black')
ax.set_title(r'r = %2.0f' %r  + '   x(0) = %2.5f' %x0, fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
#ax.set_ylim([-1, 1])
#ax.set_yticks(np.arange(y1,y2,dy))

xP = t; yP = xS; 
ax.plot(xP,yP,'b',lw = 2)    
       
fig.savefig('a2.png')





#%% CELL 3

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,3)


r1 = linspace(-1/4,1,999)
r2 = linspace(-1/4,0,999)

r1 = 0.5
z1 = 1 + (1 + 4*r1)**0.5
z2 = 1 - (1 + 4*r2)**0.5
 
xe1 = (0.5*z1)**0.5
xe2 = -xe1
xe3 = (0.5*z2)**0.5
xe4 = -xe3


y1 = 0
y2 = 4


dy = 1



fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.22,\
                    right = 0.92, hspace = 0.20,wspace=0.2)
  
ax.set_xlabel('r',color = 'black',fontsize = 12)
ax.set_ylabel('$x_e$',color= 'black',fontsize = 14)

ax.xaxis.grid()
ax.yaxis.grid()
#ax.set_ylim([y1, y2])
#ax.set_xlim([-20,1])
ax.set_xticks(np.arange(-0.5,1.1,0.25))

xP = r1; yP = xe1; col = [0,0,1]
ax.plot(xP,yP,color = col,lw = 2) 

xP = r1; yP = xe2; col = [0,0,1]
ax.plot(xP,yP,color = col,lw = 2) 

xP = r2; yP = xe3; col = [1,0,0]
ax.plot(xP,yP,color = col,lw = 2) 

xP = r2; yP = xe4; col = [1,0,0]
ax.plot(xP,yP,color = col,lw = 2) 

xP = [-0.5,0]; yP = [0,0]; col = [0,0,1]
ax.plot(xP,yP,color = col,lw = 2) 
xP = [0,1]; yP = [0,0]; col = [1,0,0]
ax.plot(xP,yP,color = col,lw = 2) 

#ax.legend() 
      
fig.savefig('a3.png')