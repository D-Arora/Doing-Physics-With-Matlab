# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 07:37:42 2023

@author: Owner
"""

# animate_D.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from matplotlib.animation import FuncAnimation


# create a figure and axes
fig = plt.figure(figsize=(12,5))
ax1 = plt.subplot(1,1,1)   

# set up the subplots as needed
ax1.set_xlim(( 0, 500))            
ax1.set_ylim((-11, 11))
ax1.set_xlabel('t')
ax1.set_ylabel('y')

# Create objects that will change in the animation. 
# Initially empty, have new values in the animation.

txt_title = ax1.set_title('')
line1, = ax1.plot([], [], 'blue', lw=2) 

numF = 20
A = np.linspace(0,1,numF)
T = 100

def drawframe(n):
    x = np.linspace(0, 500, 1000)
    y = A[n]*np.sin(2 * np.pi * x/T)
   
    line1.set_data(x, y)
    
    txt_title.set_text('Frame = {0:4d}'.format(n))
    return (line1)

# blit=True updates only the changed parts
anim = animation.FuncAnimation(fig, drawframe, frames=numF, interval=20, blit=True )



