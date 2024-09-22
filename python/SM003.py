# -*- coding: utf-8 -*-
"""
SM003.py    Sept 2024

QUANTUM MECHANICS
STATISTICAL MECHANICS  
Mawell-Boltzmann distribution: colliding disks

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH MATLAB 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSMD02A.htm
    
"""

#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np
import math
from math import factorial 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag
from scipy.integrate import odeint, quad, dblquad, simps
import pylab

tStart = time.time()

#%%  FUNCTIONS

def motion(x,y,vx,vy):
    if x < R: vx = -vx
    if x > 10-R: vx = -vx
    if y < R: vy = -vy
    if y > 10-R: vy = -vy
    x = x + vx*dt; y = y + vy*dt
    return x,y,vx,vy,

def coll(p,q):
    r1 = np.array([x[p],y[p]]);   r2 = np.array([x[q],y[q]])
    v1 = np.array([vx[p],vy[p]]); v2 = np.array([vx[q],vy[q]])
    d = np.linalg.norm(r1 - r2)
    d2 = d**2
    x1 = x[p]; x2 = x[q]; y1 = y[p]; y2 = y[q]
    if d < 1*D:
       # V1 = zeros(2); V2 = zeros(2)
        z1 = v1 - 2*m2 / M * np.dot(v1-v2, r1-r2) / d2 * (r1 - r2)
        z2 = v2 - 2*m1 / M * np.dot(v2-v1, r2-r1) / d2 * (r2 - r1) 
        v1 = z1; v2 = z2
        while d < D:
              x1 = x1 + z1[0]*dt; y1 = y1 + z2[1]*dt
              d = sqrt((x1-x2)**2 + (y1-y2)**2)
              if x1 < R: z1[0] = -z1[0]
              if x1 > 10-R: z1[0] = -z1[0]
              if y1 < R: z1[1] = -z1[1]
              if y1 > 10-R: z1[1] = -z1[1]
    return v1,v2,x1,x2


#%%  SETUP
N = 4   # number of hard disks
nF = factorial(N-1)

nT = 800  # number of time steps (number of frames for animation)

x = 1+7*np.random.random(N)    # random initial positions
y = 1+7*np.random.random(N)

vx = -6+6*np.random.random(N)   # random initial velocites
vy = -5+5*np.random.random(N)

R = 1   # radis of disk
D = 2*R   # diameter of disks
dt = 0.01  # time step
t = 0      # initial time
m1 = 1; m2 = 1; M = m1 + m2   # disk masses

# Position and velocity vectors
r1 = zeros(2); r2 = zeros(2)
v1 = zeros(2); v2 = zeros(2)

# Assign color to a disk
col = zeros([N,3])
col[0,:] = [1,0,0]; col[1,:] = [0,0,1]
col[2,:] = [1,0,1]; col[3,:] = [0,0,0]

# Vector for position of 4 disks
G = [0,0,0,0]


#%% SETUP plot for animation
plt.rcParams["figure.figsize"] = (4,4)
fig1, axes = plt.subplots()
axes.set_aspect('equal', adjustable='box')
axes.set_aspect(1)
axes.set_xlim([0,10])
axes.set_ylim([0,10])
axes.xaxis.set_ticks([])
axes.yaxis.set_ticks([])

for c in range(4):
    G[c], = axes.plot([], [], 'o', color = col[c,:],ms = 30) 

def init():
   for c in range(4): 
       G[c].set_data([], [])
   return  G 


#%% ANIMATION CODE
def animate(n):  
    
# WALL    
         for k in range(4):
            x[k],y[k],vx[k],vy[k] = motion(x[k],y[k],vx[k],vy[k])       
    
# COLLISIONS
            p = 0; q = 1
            for c in range(nF):
                v1, v2, x1,x2 = coll(p,q)
                x[p] = x1; x[q] = x2        
                vx[p] = v1[0]; vy[p] = v1[1]
                vx[q] = v2[0]; vy[q] = v2[1] 
                q = q+1
                if q > N-1: p = p+1; q = p+1
            
# WALL    
            for k in range(4):
             x[k],y[k],vx[k],vy[k] = motion(x[k],y[k],vx[k],vy[k])            
         
         
         for c in range(4):
             x1 = x[c]; y1 = y[c]
             G[c].set_data([x1], [y1])
                        
         time.sleep(0.001)
         
         return  G #G[0], G[1],G[2],G[3]

anim = FuncAnimation(fig1, animate, init_func = init, 
                     frames = nT, interval = 20, blit = True, repeat = False)

# anim.save('AAA.gif', fps = 10)   
    

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)