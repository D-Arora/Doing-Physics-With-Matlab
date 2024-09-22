# -*- coding: utf-8 -*-
"""
qmSM02.py    Aug 2024

QUANTUM MECHANICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSM01.pdf


STATISTICAL MECHANICS: MAXWELL - BOLTZMANN DISTRIBUTION
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
import random
import sys

tStart = time.time()


#%%
num = 5000
E = 9
box = zeros([num,E]); 
Esum = 0
B = np.arange(0,9,1)

for k in range(num):
    BOX = zeros(E); S = 0
    for c in range(1500):
        s = random.randint(0,E-1)
        S = S + s
        if S < 9:
           BOX[s] = BOX[s] + 1
        if S > 8:
           S = S - s   
        if S == 8:
           c = 100
    box[k,:] = BOX
    print(box[k,:])
    
N = np.sum(box,axis = 0)     
print(N)
  


#%% Figure 1: Plot of Maxwell speed distribution
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

# ax.set_xlabel('v [ m/s ]',fontsize = 12)
# ax.set_ylabel('f$_M$  [ x10$^{-3}$  s/m ]',fontsize = 12)
# ax.xaxis.grid()
ax.set_ylim([0,10000])
#xP = v
# yP = fM*1e3
for c in range(9):
    xP = [c,c]
    yP = [0,N[c]]
    ax.plot(xP,yP,'b',lw = 3)




#%% SAVE FIGURES
#  fig1.savefig('a1.png')
#  fig2.savefig('a2.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)



