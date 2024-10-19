# -*- coding: utf-8 -*-

"""
op004.py    oct 2024

COMPUTATIONAL OPTICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/op004.pdf

COMPUTATION OF DOUBLE (SURFACE) INTEGRALS

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag, real
import random
from pylab import imshow,show,gray,jet,colorbar,xlabel,ylabel,title,hsv,hot
from scipy.integrate import odeint, quad, dblquad, simps

tStart = time.time()


#%% INPUTS >>>
# Integration limits
x1 = -1; x2 = 1
y1 = 0; y2 = 1

def Fn(xx,yy):
    F = xx**2 + xx**2
    return F
    
#%% METHOD 1: Integration a set of [1D] strips byf Simpson's [1D] rule 
# XY grid points and [2D] grid 
Nx = 111
Ny = 111

x = linspace(x1,x2,Nx+1)
y = linspace(y1,y2,Ny+1)
dy = y[2]-y[1]
X,Y = np.meshgrid(x,y)
# Function to be integrated
F = Fn(X,Y)
# Initialize irradiance
I1 = 0
# Compute [2D] integral
for c in range(Ny):
    fn = F[c,:]
    I1 = I1 + simps(fn,x)*dy

print('Method 1')
print('irradiance I = %0.6f' %I1)    


#%% METHOD 2: [2D] Simpson's rule
# [2D] grid
num = 55
x1 = -1
x2 = 1
y1 =0
y2 = 1

x = linspace(x1,x2,num)
y = linspace(y1,y2,num)
hx = (x2-x1)/(num-1); hy = (y2-y1)/(num-1)
h = hx * hy / 9
X,Y = np.meshgrid(x,y)
# Function to be integrated
F = Fn(X,Y)
# Compute [2D] integral: [2D] Simpson coefficients and integral
sc = 2*ones(num)
sc[1::2]= 4
sc[0] = 1
sc[num-1] = 1
sc1,sc2 = np.meshgrid(sc,sc)
S = sc1*sc2
SF = S*F
# Compute [2D] integral
I2 = h * np.sum(SF)

print('  ')
print('Method 2')
print('irradiance I = %0.6f' %I2)    


#%%
# # # https://docs.scipy.org/doc/scipy/tutorial/integrate.html

# # # https://sites.millersville.edu/rbuchanan/math375/DoubleInt.pdf

# from scipy import integrate
# def f(x, y):
#     return x**2 + y**2

# def bounds_y():
#     return [0, 1]

# def bounds_x(y):
#     return [-1, 1]

# II = integrate.nquad(f, [bounds_x, bounds_y])
# print(II)

