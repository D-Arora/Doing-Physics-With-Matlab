# -*- coding: utf-8 -*-

# turningPoint.py


# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
import math
from numpy import pi, sin, cos

t = np.linspace(0,10,199)
y = sin(2*pi*t)


# Findpeaks
L = len(y)
a1 = y[0]; a2 = y[1]
flag = 2
indexMax = np.zeros(20)

if a2 > a1:
    flag = 1
v = 0
for c in range(L-1):
    y1 = y[c]     
    y2 = y[c+1]

    if flag == 1 and y2 > y1:
       c = c+1
    if (flag == 1 and y2 < y1):
       indexMax[v] = c
       v = v + 1
       c = c+1
    if y2 <= y1:
      flag = 0
    if y2 > y1:
      flag = 1

ind = indexMax[indexMax > 0]
ind = ind.astype(int)


font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

# Fig 1  t vs theta  ---------------------------------------------------------
fig = plt.figure(figsize = (4, 3))

fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                     right = 0.96, hspace = 0.2,wspace=0.2)

xP = t; yP = y
plt.plot(xP,yP,linewidth=2,color='r')
xP = t[ind]; yP = y[ind]
plt.plot(xP,yP,'bo')
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
plt.grid('visible')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'$y$', fontdict = font1)

