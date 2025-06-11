# -*- coding: utf-8 -*-
"""
S004A.py               06 June 2025
 
 RAY OPTICS
   Matrix methods in paraxial optics
   TRANSFOMATION MATRICES: ABCD matrices
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/S003.pdf

Functions to generate the transformation matrices ABCD

THIN  LENS
      1/s0 + 1/s1 = 1/f  m = -s1/s0 --> 
      plots of s1 and m against s0
"""

import numpy as np
from numpy import arange, exp, linspace, zeros, amax, sqrt
from numpy import pi, sin, cos, tan, arctan
import matplotlib.pyplot as plt
import time

tStart = time.time()

plt.close('all')


#%% THIN LENS EQUATION
# Focal length
f = 5
# Object plane 
sA = 0.1; sB = 20; N = 999
s0 = linspace(sA,sB,N)
# Image plane
s1 = s0*f/(s0-f)
# Magnification
mag = -s1/s0

plt.rcParams["figure.figsize"] = (5,3)
plt.rcParams['font.size'] = 12
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax1 = ax.twinx()
X = -f; Y = 0; ax.plot(X,Y,'ko',ms = 6)
X = -s0; Y = s1; ax.plot(X,Y,'b',lw = 2)   
X = -2*f; Y = 2*f; ax.plot(X,Y,'bo',ms = 5)
ax.set_xlabel('s$_0$', fontsize=14)
ax.set_ylabel('s$_1$', fontsize=14, color = 'b')
ax.set_title('f = %0.2f' %f, fontsize = 12)
ax.tick_params(axis='y', colors='blue')
ax.set_ylim([-40,40])
ax.grid()

X = -s0; Y = mag;ax1.plot(X,Y,'r',lw = 1) 
X = [-20,0]; Y = [-1,-1]; ax1.plot(X,Y,'r',lw = 1) 
X = -2*f; Y = -1; ax1.plot(X,Y,'ro',ms = 5)

ax1.set_ylim([-5,5])
ax1.set_ylabel('mag', fontsize=12, color = 'r')
ax1.tick_params(axis='y', colors='red')

fig2.tight_layout()



#%% Save figures
# fig2.savefig('a2.png')   
      
      


