# -*- coding: utf-8 -*-
"""
Created on Sun Aug 10 10:16:11 2025

@author: Owner
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace

plt.close('all')

tStart = time.time()


#%%  SOLVE ODE x
def lorenz(t, state):    
    x = state
    dx = sin(x)
    return dx 

x0 = pi/2
tMax = 10

5
t = linspace(0,tMax,9999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 

plt.plot(t,xS)