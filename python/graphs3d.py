# -*- coding: utf-8 -*-

"""
qm061.py    June 2024

Ian Cooper 
        matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm061.pdf

QUANTUM MECHANICS     ANGULAR MOMENTUM:  AZIMUTHAL WAVEFUNCTION
  
"""

# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, radians
import scipy.special 
import scipy.special as sps
from scipy import *
import scipy.linalg as la
import math
import cmath
import os
import matplotlib.pyplot as plt
from pylab import *
from math import *
from numpy import linspace,matrix,array
from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve
from scipy.integrate import odeint, solve_ivp, simps
from scipy.sparse.linalg import eigsh, eigs #Solves the Eigenvalue problem
from scipy.sparse import diags #Allows us to construct our matrices
from matplotlib.animation import FuncAnimation, PillowWriter 
from scipy.special import sph_harm
import time

tStart = time.time()


#%%   [3D] line plots
# Creating an empty figure
plt.rcParams["figure.figsize"] = (4,3)
plt.rcParams['font.size'] = 10

fig1 = plt.figure()
fig1.subplots_adjust(top = 0.88, bottom = 0.10, left = 0.05,
                    right = 0.90, hspace = 0.20,wspace=0.20)
# Defining the axes as 3D
ax = plt.axes(projection="3d")
# Defining x,y, z values
z = np.linspace(0, 1, 100)
x = z * np.cos(25 * z)
y = z * np.sin(25 * z)
# [3D] line plot
ax.plot3D(x, y, z,'b',lw = 2)
ax.set_title('A simple [3D] line plot')
ax.set_xlabel('x',fontsize=12)
ax.set_ylabel('y',fontsize=12)
ax.set_zlabel('z',fontsize=12)
ax.set_xticks(arange(-1,1,0.5))
ax.set_yticks(arange(-1,1,0.5))

ax.zaxis.set_major_formatter('{x:.02f}')


#%%   SURFACE PLOT
fig2 = plt.figure()
 
# Creating our empty subplots
ax = fig2.add_subplot(1, 1, 1, projection='3d')
 
# Creating our data points
x1= arange(-1,1,0.1)
y1= arange(-1,1,0.1)
 
# Creating a mesh grid
x1,y1= meshgrid(x1,y1)
 
# Creating a cosine function with the
# range of values from the meshgrid
z1= cos(x1* np.pi/2)*sin(y1* np.pi/2)
 
# Creating a wireframe plot with the points
# x1,y1,z1 along with the plot line as red
ax.plot_surface(x1, y1, z1, color="yellow")
ax.set_title('Surface plot')



#%% INPUTS
# oribtal quantum nuumber L
L = 1
# magnetic quantum number mL
mL = -1
# Grid point N
N = 99

 
#%% SETUP
# Polar angle  [rad]  / Azimuthal angle  [ra]
theta = linspace(0,pi, N)
phi   = linspace(0,2*pi, N)
# [D] meshgrid
THETA, PHI =  np.meshgrid(theta, phi) 
# Calculate the Cartesian coordinates of each point in the mesh
xyz = np.array([sin(THETA)*sin(PHI), sin(THETA)*cos(PHI), cos(THETA)])
# Calculate spherical harmonic
Y = sph_harm(abs(mL), L, PHI, THETA)
YR = Y.real; YI = Y.imag

if mL < 0:
     Y = sqrt(2) * (-1)**mL * Y.imag
elif mL > 0:
     Y = sqrt(2) * (-1)**mL * Y.real

Yx, Yy, Yz = abs(Y) * xyz


plt.rcParams["figure.figsize"] = (2,2)
fig3 = plt.figure()
fig3.subplots_adjust(top = 0.8, bottom = 0.0, left = 0.0,
                    right = 1, hspace = 0.20,wspace=0.20)
 
# Creating our empty subplots
ax = fig3.add_subplot(1, 1, 1, projection='3d')
 

#ax.plot_surface(Yx, Yy, Yz, color="blue")
ax.plot_surface(Yx, Yy, Yz, cmap='jet')
ax.set_title('l = %2.0f' %L + '   $m_l = %2.0f$' %mL)
L2 =0.2; L1 = -L2
ax.plot([L1,L2], [0,0], [0,0], c='0.4', lw=1, zorder=10)
ax.plot([L1,L2], [0,0], [0,0], c='0.4', lw=1, zorder=10)
ax.plot([0,0], [L1,L2], [0,0], c='0.4', lw=1, zorder=10)
ax.plot([0,0], [0,0], [L1,L2], c='0.4', lw=1, zorder=10)
ax.axis('off')

fig3.savefig('a1.png',bbox_inches='tight',dpi=800,
             pad_inches=0.1, transparent=True)

# Displaying the plot
#plt.show()

# #%%
# phi = linspace(0,2*pi,599)         # arimuthal angle  [rad]

# # r real part / im imaginary part / pd probability density

# def graph(R,C,mL):
#     r = abs(cos(mL*phi)); im = abs(sin(mL*phi)); pd = r**2 + im**2
#     ax[R,C].plot(phi, r,'b')
#     ax[R,C].plot(phi, im,'k')
#     ax[R,C].plot(phi, pd,'r',lw = 4)
#     ax[R,C].set_rmax(1)
#     ax[R,C].set_rticks([])  # Less radial ticks
#     ax[R,C].grid(True)
#     ax[R,C].set_xticks(arange(0,2*pi, pi/6))
#     ax[R,C].set_title('m$_L$ = %2.0f' % mL)


# #%%
# plt.rcParams['font.size'] = 10
# plt.rcParams["figure.figsize"] = (7,7.2)
# fig1, ax = plt.subplots(nrows=3, ncols=2,subplot_kw={'projection': 'polar'})
# #fig1, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# fig1.subplots_adjust(top = 0.95, bottom = 0.07, left = 0.120,\
#                     right = 0.98, hspace = 0.60,wspace=0.40)

# R = 0; C = 0; mL = 0 
# graph(R,C,mL)
# R = 0; C = 1; mL = 1  
# graph(R,C,mL)
# R = 1; C = 0; mL = 2  
# graph(R,C,mL)
# R = 1; C = 1; mL = 3  
# graph(R,C,mL)
# R = 2; C = 0; mL = 4 
# graph(R,C,mL)
# R = 2; C = 1; mL = 5  
# graph(R,C,mL)

# #fig1.savefig('a1.png')




# #%%
# tExe = time.time() - tStart
# print('  ')
# print('Execution time = %2.0f s' % tExe)


# #%%
# from scipy.special import clpmn, lpmn
# z = pi
# L = 1
# mL = 0
# zz = lpmn(mL,L,z)


# A = list(zip(*zz))

# res = np.array(A, dtype=float)

# B = res[0,0,0]
# print(B)


#https://colab.research.google.com/github/caiociardelli/sphglltools/blob/main/doc/3_Associated_Legendre_polynomials.ipynb


