# -*- coding: utf-8 -*-

# op003.py               March 2024
 
#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/op001.pdf

# [1D] integration

import numpy as np
from numpy import pi, sin, cos, exp, linspace 
from numpy.linalg import eig
from scipy.integrate import odeint, quad, dblquad, simps
from scipy import pi, sqrt
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import axes3d
import sympy as sym


def func(x):
    #f = cos(x)
    f = exp(2+x*(1j))
    return f

# >>>>> Inputs: Grid points,limits
N = 299
a = 0; b = pi/2

# function f(x)
x =linspace(a,b,N)
f = func(x)

# Integrate function uisng quad
Iquad, Ierr = quad(func,a,b)
print(Iquad,Ierr)

# Integrate function using simps
Is = simps(f,x)
print(Is)

# Plot function
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)

# Figure 1    t vs x
fig1, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('real(f)',color = 'black',fontsize = 12)
axes.xaxis.grid()
axes.yaxis.grid()
xP = x; yP = f
axes.plot(xP, yP,'b',lw = 2)
fig1.tight_layout()
# fig1.savefig('a1.png')


#%%  Symoblic integration
x = sym.Symbol('x')
y = sym.Symbol('y')

y = 6 * x ** 5
y = sym.sin(x)
y = sym.log(x)
y = 2 * x + sym.sinh(x)
sym.integrate(y, x)

# Highlight the code and use F9 to execute the code
# It is possible to compute definite integral:
# y = x**4;  a = -1; b = 1; sym.integrate(y, (x, a, b))

# sym.integrate(sym.sin(x), (x, 0, sym.pi))

# sym.integrate(sym.cos(x), (x, -sym.pi / 2, sym.pi / 2))

# Also improper integrals are supported as well:
# sym.integrate(sym.exp(-x), (x, 0, sym.oo))

# sym.integrate(sym.exp(-x ** 2), (x, -sym.oo, sym.oo))

# https://scipy-lectures.org/packages/sympy.html#integration




