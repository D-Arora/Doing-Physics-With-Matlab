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
    xDot = r + x**2
   
    return xDot

#%%  SOLVE ODE
# Solve ODE for x
def lorenz(t, state):    
    xS = state
    dx = r + xS**2
    return dx  

def graph(xP,yP,col):
    ax.plot(xP,yP,color = col,lw = 2)
    return
#%%  CELL 1   
r1 = -16; r2 = 0; r3 = 16
x = linspace(-8,6,299)

xDot1 = funct(r1,x)
xDot2 = funct(r2,x)
xDot3 = funct(r3,x)

xe1 = [-(-r1)**0.5, (-r1)**0.5] 
xe2 = [-(-r2)**0.5, (-r2)**0.5] 
xe3 = [-(-r3)**0.5, (-r3)**0.5]

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)

fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.18,\
                    right = 0.95, hspace = 0.20,wspace=0.2)

ax.set_xlabel('x',color = 'black')
ax.set_ylabel('xDot',color= 'black') 
   
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_xlim([-8, 8])
ax.set_xticks(np.arange(-8,8.1,2))
ax.set_ylim([-20, 20])
ax.set_yticks(np.arange(-20,21,5))
ax.set_title(r'r = %2.0f' %r1,  fontsize = 12)

# r < 0
# xP = x; yP = xDot1; col = [0,0,1]
# ax.plot(xP,yP,color = col,lw = 2, label = 'r = -16')
# xP = xe1[1]; yP = 0; col = [1,0,0]
# ax.plot(xP,yP,'o',color = col, ms = 8)
# xP = xe1[0]; yP = 0; col = [0,0,1]
# ax.plot(xP,yP,'o',color = col, ms = 8)

# ax.arrow(-7,0,1,0,lw = 2,color = 'b', width = 0.2, head_width = 3,
#          head_length = 0.4)

# ax.arrow(3,0,-1,0,lw = 2,color = 'b', width = 0.2, head_width = 3,
#          head_length = 0.4)

# ax.arrow(5,0,1,0,lw = 2,color = 'b', width = 0.2, head_width = 3,
#          head_length = 0.4)


r = 0 
ax.set_title(r'r = %2.0f' %r,  fontsize = 12)
ax.set_xlim([-6, 6])
ax.set_xticks(np.arange(-6,6.1,2))
ax.set_ylim([-5, 20])
ax.set_yticks(np.arange(-5,21,5))
xP = x; yP = xDot2; col = [0,0,1]
ax.plot(xP,yP,color = col,lw = 2, label = 'r = 0')

xP = xe2[0]+0.1; yP = 0; col = [1,0,0]
ax.plot(xP,yP,'o',color = col, ms = 8)
xP = xe2[0]-0.1; yP = 0; col = [0,0,1]
ax.plot(xP,yP,'o',color = col, ms = 8)

ax.arrow(2,0,1.2,0,lw = 2,color = 'b', width = 0.2, head_width = 2,
          head_length = 0.8)

ax.arrow(-4,0,1.2,0,lw = 2,color = 'b', width = 0.2, head_width = 2,
          head_length = 0.8)

# r > 0
#xP = x; yP = xDot3; col = [1,0,1]
#ax.plot(xP,yP,color = col,lw = 2, label = 'r = +16')






# ax.legend()

fig.savefig('a1.png')


#%% CELL 2
r = 0
x0 = -8
y1 = -10
y2 = 2
dy = 2

N = 9999
t1 = 0
t2 = 20
t = linspace(t1,t2,N)

sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 

#%%

fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.90, bottom = 0.23, left = 0.18,\
                    right = 0.98, hspace = 0.20,wspace=0.2)
  
ax.set_xlabel('t',color = 'black')
ax.set_ylabel('x',color= 'black')
ax.set_title(r'r = %2.0f' %r  + '   x(0) = %2.5f' %x0, fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([y1, y2])
ax.set_yticks(np.arange(y1,y2,dy))

xP = t; yP = xS; col = [0,0,1]
ax.plot(xP,yP,'b',lw = 2)    
       
fig.savefig('a2.png')


#%% CELL 3

R = linspace(-20,0,599)
xe1 = (-R)**0.5
xe2 = -xe1
y1 = -5
y2 = 5
dy = 1

fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.25,\
                    right = 0.92, hspace = 0.20,wspace=0.2)
  
ax.set_xlabel('r',color = 'black',fontsize = 12)
ax.set_ylabel('$x_e$',color= 'black',fontsize = 14)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([y1, y2])
ax.set_xlim([-20,1])
ax.set_xticks(np.arange(-20,0.1,4))

xP = R; yP = xe1; col = [1,0,0]
ax.plot(xP,yP,color = col,lw = 2,label = 'unstable') 
xP = R; yP = xe2; col = [0,0,1]
ax.plot(xP,yP,color = col,lw = 2,label = 'stable') 
xP = -16; yP = -4; col = [0,0,1]
ax.plot(xP,yP,'o',color = col,ms = 8) 
xP = -16; yP = 4; col = [1,0,0]
ax.plot(xP,yP,'o',color = col,ms = 8) 
ax.legend(fontsize = 10) 
       
fig.savefig('a3.png')