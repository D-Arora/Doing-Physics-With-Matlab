# -*- coding: utf-8 -*-

# op002.py          March 2024

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/op001.pdf

# [2D] integration

# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import axes3d

tStart = time.time()

# >>>>> Input number of grid points: N must be odd
N = 999

# p rho /  q phi

# >>>>> Setup RP meshgrid
p1 = 0;  p2 = 2; q1 = 0; q2 = 2*pi
p = linspace(p1,p2,N)
q = linspace(q1,q2,N)
pp, qq = np.meshgrid(p,q)

# Mask matrix M
M = np.ones([N,N])

# Function 
f = pp


# Simpson [2D] coefficients
S = np.ones(N)
R = np.arange(1,N,2);   S[R] = 4;
R = np.arange(2,N-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy

# Calculate integral
hp = (p2-p1)/(N-1); hq = (q2-q1)/(N-1)
h = hp * hq / 9
integral = h*sum(sum(pp*f*S))

print('\n Integral  =  ',integral)


xx = pp*cos(qq); yy = pp*sin(qq)

#%% GRAPHICS
plt.rcParams['font.size'] = 10
fig = plt.figure(figsize=(4,3))
ax = plt.axes(projection='3d')
ax.plot_surface(xx, yy, f, cmap='jet',
      edgecolor='none', alpha=1,antialiased=True)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.set_zlabel('f', fontsize=12)
fig.tight_layout()
ax.view_init(17,-110,0)
# fig.savefig('a1.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

