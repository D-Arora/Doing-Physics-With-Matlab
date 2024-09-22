# -*- coding: utf-8 -*-
"""
SM003.py    Sept 2024

QUANTUM MECHANICS
STATISTICAL MECHANICS  
Mawell-Boltzmann distribution: 
    Animation: molecular motion of gas molecules

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH MATLAB 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSMS02A.htm
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


#%%  ANIMATION SETUP
plt.rcParams["figure.figsize"] = (4,4)
fig1, axes = plt.subplots()
axes.set_aspect('equal', adjustable='box')
axes.set_aspect(1)
axes.set_xlim([0,10])
axes.set_ylim([0,10])
axes.xaxis.set_ticks([])
axes.yaxis.set_ticks([])

G, = axes.plot([], [], 'bo',ms = 3) 

def init():
       G.set_data([], [])
       return  G, 


#%% ANIMATION CODE
def animate(n):  
    x, y = 10*np.random.random([2,1000]) 
    G.set_data([x], [y])
    time.sleep(0.1)
    return  G, 

anim = FuncAnimation(fig1, animate, init_func = init, 
                     frames = 80, interval = 5, blit = True, repeat = False)

#  anim.save('AAA.gif', fps = 10)   
    

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)