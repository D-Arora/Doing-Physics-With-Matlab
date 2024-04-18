# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 20:17:35 2024

@author: Owner
"""

import time
import math
import numpy as np
import pylab as py
from matplotlib import pyplot as plt 
from numpy import linspace,sin,cos,exp, zeros, ones, pi, diag

#%% FUNCTIONS

def firstDer(N,dx):
    v  = ones(N-1)
    M1 = diag(-v,-1)
    M2 = diag(v,1)
    M = M1+M2
    M[0,0] = -2; M[0,1] = 2; M[N-1,N-2] = -2; M[N-1,N-1] = 2
    MF = M/(2*dx) 
    return MF

def secondDer(N,dx):
     v = -2*ones(N)
     M1 = np.diag(v)
     v = np.ones(N-1)
     M2 = np.diag(v,1)
     M3 = np.diag(v,-1)
     M = M1+M2+M3
     M[0,0] = 1; M[0,1] = -2; M[0,2] = 1
     M[N-1,N-3] = 1; M[N-1,N-2] = -2; M[N-1,N-1]=1
     MS = M/(dx**2) 
     return MS


#%% Function to be differentiated
N = 199                  # Number of array elements
x1 = 0; x2 = 10          # Function limits
x = linspace(x1,x2,N)    # x array
dx = x[2]-x[1]

y  = sin(2*x)              #  Function to be differentiated 

MF = firstDer(N,dx)      #  1st derivative matrix

FD = MF@y                # 1st derivative

MS = secondDer(N,dx)     # 2nd derivative matrix

SD = MS@y                # Second derivative


#%%  GRAPHICS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,4)

fig, ax = plt.subplots(1)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('y, dy/dx, -d$^2$y/dx$^2$ ',color= 'black')
ax.set_xlabel('x ',color = 'black')
ax.set_xlim([0, 10])
#ax.set_ylim([-1.1, 1.1])
#ax.set_yticks(np.arange(-1,1.2,0.50))
fig.tight_layout()
ax.plot(x,y,'b', lw = 3, label = 'y')
ax.plot(x,FD,'r',lw = 2, label = 'dy/dx')
ax.plot(x,-SD,'k',lw = 1, label = '-d$^2$y/dx$^2$ ' )
ax.legend()

# fig.savefig('a1.png')




