
# https://www.oreilly.com/library/view/python-data-science/9781491912126/ch04.html

import numpy
import numpy as np

import scipy.special 
import scipy.special as sps
import math
import os
import matplotlib.pyplot as plt
pi = np.pi
from pylab import figure, plot,show,legend,xlabel,ylabel,title,grid
from pylab import ylim,xlim
from numpy import linspace,sin,cos,exp


x1 = 0; x2 = 100; Nx = 299;
L = 20; A1 = 50; A2 = 25
x = linspace(x1,x2,Nx)
y1 = A1*sin(2*pi*x/L)
y2 = A2*cos(2*pi*x/L)

plt.rcParams["figure.figsize"] = (5,4)

#fig = plt.figure(figsize = (5, 4))
fig, ax = plt.subplots(2)

#ax[0] = fig.add_axes([0.1, 0.5, 0.8, 0.4])
#ax[0].set_figwidth(4)
#ax.figuresize(4,4)
#ax[0] = fig.add_axes([0.1, 0.5, 0.8, 0.4], \
#                          xticklabels=[], ylim=([-50, 50])
#ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.4],\
#                          ylim=(-1.2, 1.2))

ax[0].plot(x, y1,'r',label = 'sine')
ax[0].plot(x, y2,'b',label = 'cos')
ax[0].set(xlim =[0, 100], yticks=[-50,50])
ax[0].set(xticks=[0,100], yticks=[])
#ax[0].xlim([0,100])
ax[0].axis('equal')
ax[0].legend(loc ='upper right',frameon = True, ncol = 2)

ax[1].plot(x, y2)




