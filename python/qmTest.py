# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 16:06:33 2024

@author: Owner
"""

# https://phys.iit.edu/~segre/phys405/12F/lecture_03.pdf


# qmTest.py

import time
import math
import numpy as np
import pylab as py
from matplotlib import pyplot as plt 
from numpy import linspace,sin,cos,exp, zeros, ones, pi, diag, sqrt
from scipy.integrate import odeint, quad, dblquad, simps


N = 5

v = np.arange(0,N,1)

M1 = diag(v)

y = v

M2 = M1@y

M3 = y@M1@y

print(M2)

print(M3)

#%%

x = linspace(0,1,9)
y = sqrt(2)*sin(pi*x)


X = diag(x)
XAVG = y@X@y
xavg = simps(y*y,x)
print(XAVG)

print(xavg)

M1 = y@y


print(M1)