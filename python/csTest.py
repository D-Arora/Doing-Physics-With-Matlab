# -*- coding: utf-8 -*-
"""
csTest
"""

# LIBRARIES  ================================================================
import numpy as np
from numpy import pi, sin, cos, linspace
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')

#%%
N = 15
t = linspace(0,10,N)
x = linspace(0,2,N)

f = x*(1-x)


T,X = np.meshgrid(t,x)

dX = np.ones([N,N])
F = X*(1-X)
dY = F/(np.sqrt(dX**2 + F**2))
dX = dX/(np.sqrt(dX**2 + F**2))

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_ylim([0,2])
ax.set_xlim([0,10])

ax.quiver(T,X,dX,dY,color = 'k')


fig4.tight_layout()  

