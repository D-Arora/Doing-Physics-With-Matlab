# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 15:07:29 2022

@author: Owner
"""

import numpy
import numpy as np

# Special functions
# https://docs.scipy.org/doc/scipy/reference/special.html


import scipy.special 
import scipy.special as sps

import math

import os
import matplotlib.pyplot as plt

# =====================================================================

pi = np.pi
from pylab import figure, plot,show,legend,xlabel,ylabel,title,grid
from pylab import ylim,xlim
from numpy import linspace,sin,cos,exp

# =====================================================================
T = 3
wL = 20
x = linspace(0,100,201)
t = linspace(0.0, 6, 201)
y = sin(2*pi*x/wL - 2*pi*t[1])/T

#figure(1)
plt.subplots(figsize=(10, 5))
plt.rcParams['font.size'] = '14'

for c in range(201):
    plt.clf()
    z = 0.8*sin(2*pi*x/wL - 2*pi*t[c]/T)
    plot(x,z,linewidth = 3,color = 'b')
    xlabel('x',fontsize = 16)
    ylabel('y',fontsize = 16)
    title('t = {:2.2f}'.format(t[c]))
    ylim(-1,1)
    xlim(0,100)
    grid(color='k', linestyle='-', linewidth = 0.2)
    plt.pause(0.001)
#    print('pi = {:2.5f}     e = {:2.4f}   sum = {:3.2e}  '.format(pi,e, S))
    
