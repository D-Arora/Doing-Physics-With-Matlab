# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 14:00:38 2024

@author: Owner
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
epsi = 0.01
# Declare the model
f_dir = np.arange(0,1.3,0.001)
A_s = np.zeros(len(f_dir))

def myModel(y, t, f):
    dy0 = y[1]
    dy1 = -epsi*y[1]-np.sin(y[0]) - f*np.cos((1.01)*t)*np.cos(y[0])
    return [dy0, dy1]

i = 0
for f in f_dir:
    time = np.arange(0.0, 2000,0.01)
    yinit = np.array([np.pi/2, 0])
    y = odeint(myModel, yinit, time, args=(f,))
    A_s[i] = np.abs(np.max(y[-600:-1,0])- np.min(y[-600:-1,0]))
    i += 1


plt.plot(f_dir,A_s,',')
plt.xlabel(r'$f_s$')
plt.ylabel(r'$A_s$')

plt.show()