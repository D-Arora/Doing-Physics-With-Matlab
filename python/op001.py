# -*- coding: utf-8 -*-

# op001.py             March 2024

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
N = 2999

# >>>>> Setup xy meshgrid
#ax = 0; bx = 2; ay = 1; by = 5
#ax = 0; bx = 2; ay = 1; by = 5
#ax = 0; bx = 2; ay = 0; by = 2*pi
ax = -1; bx = 1; ay = -1; by = 1; a = 1
x = linspace(ax,bx,N)
y = linspace(ay,by,N)
xx, yy = np.meshgrid(x,y)

# Mask matrix M
M = np.ones([N,N])
#M[yy >= 1 - xx] = 0
M[xx**2+yy**2>1]=0
#M[x**2+y**2>=a**2] = 0

# Function 
#f = xx**2*yy**3
#f = 3*xx**1*yy**0
f = cos(xx)*sin(yy)

f = np.real(a**2 - xx**2 - yy**2)
#f[(xx**2 + yy**2) >= a**2] = 0
#f = f**0.5
#f = 6*np.ones([N,N])
#f = (xx**2+yy**2)
f = f*M

# Simpson [2D] coefficients
S = np.ones(N)
R = np.arange(1,N,2);   S[R] = 4;
R = np.arange(2,N-1,2); S[R] = 2
scx, scy = np.meshgrid(S,S)
S = scx*scy

# Calculate integral
hx = (bx-ax)/(N-1); hy = (by-ay)/(N-1)
h = hx * hy / 9
integral = h*sum(sum(f*S))

print('\n Integral  =  ',integral)



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
ax.view_init(14,-100,0)
# fig.savefig('a1.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

