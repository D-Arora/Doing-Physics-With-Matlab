# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 07:53:49 2023

@author: Owner
"""

# animate_E.py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import time
#from matplotlib.animation import FuncAnimation


# create a figure and axes
fig = plt.figure(figsize=(3,3))
fig.subplots_adjust(top=0.85, bottom = 0.24, left = 0.28,\
                    right = 0.90, hspace = 0.60,wspace=0.5)
ax1 = plt.subplot(1,1,1)   
#ax2 = plt.subplot(1,2,2)
# set up the subplots as needed
ax1.set_xlim(( 0, 2))            
ax1.set_ylim((-11, 11))
ax1.set_xlabel('Time')
ax1.set_ylabel('Magnitude')

txt_title = ax1.set_title('')
line1, = ax1.plot([], [], 'blue', lw=3)     
# ax.plot returns a list of 2D line objects
line2, = ax1.plot([], [], 'tan', lw=2)
#pt1, = ax2.plot([], [], 'k.', ms=20)
#line3, = ax2.plot([], [], 'grey', lw=2)
#ax1.legend(['sin','cos']);
plt.grid('visible')
numF = 100
A = np.linspace(0,10,numF)

def drawframe(n):
    x = np.linspace(0, 2, 1000)
    y1 = A[n]*np.sin(2 * np.pi * x)
   # y2 = np.cos(2 * np.pi * (x - 0.01 * n))
    line1.set_data(x, y1)
    #line2.set_data(x, y2)
    #line3.set_data(y1[0:50],y2[0:50])
    #pt1.set_data(y1[0],y2[0])
   # txt_title.set_text('Frame = {0:4d}'.format(n))
    txt_title.set_text('A = {:.2f}'.format(A[n]))
    return (line1,line2)

    #return (line1)
# blit=True updates only the changed parts
anim = animation.FuncAnimation(fig, drawframe, frames=numF, interval=50, blit=False, \
                               repeat = False)

#anim.save("ag_E.gif", dpi=250, fps=20)
#anim.save('ag_C.gif', fps=10)