# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 15:39:37 2023

@author: Owner
"""

# animate_A.py
#  ANIMATIONS
# 231226

# https://pythonguides.com/matplotlib-update-plot-in-loop/



# LIBRARIES  ================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from matplotlib.animation import FuncAnimation

T = 100;
t = np.linspace(0,500,299)

A = np.linspace(0,10,20)

y = A[2]*np.sin(2*np.pi*t/T)


font1 = {'family':'Times New Roman','color':'black','size':14}
plt.rcParams['font.family'] = ['Times New Roman']
plt.rcParams['font.size'] = 14





#x = t

figure, ax = plt.subplots(figsize=(4,3))  
plt.subplots_adjust(top=0.90, bottom = 0.20, left = 0.20,\
                     right = 0.96, hspace = 0.60,wspace=0.5)

# GUI
plt.ion()
#  Plot
plot_1, = ax.plot(t, y,'b', lw  = 2)   

xP = t
global k
k = 10
def animate_plot (c):
    global k
    for c in range(20):
        
     k = k + c   
     xP = t   
     yP = A[k]*np.sin(2*np.pi*t/T)
     plot_1.set_xdata(xP)
     plot_1.set_ydata(yP)
    
     figure.canvas.draw()
     figure.canvas.flush_events()
    #plt.plot(xP,yP,linewidth=2,color='b')    
     plt.grid('visible')
     plt.ylabel('y', fontdict = font1)
     plt.xlabel('t ', fontdict = font1)   
     plt.ylim([-11,11])
     plt.grid('visible')
     plt.draw()
     plt.pause(1) 
     return plot_1,
 
# Animated Function

ani = FuncAnimation(figure,
                    animate_plot,
                    frames=33,
                    interval=100)

# Save as gif

#ani.save('ag_001.gif', fps=10)   
    
    
# Create subplots



# Data Coordinates
# x = np.linspace(0, 20, 80)
# y = np.sin(x)


# # Labels
# plt.xlabel("X-Axis",fontsize=18)
# plt.ylabel("Y-Axis",fontsize=18)

# for value in range(150):
#     update_y_value = np.sin(x-2.5*value)
    
#     plot1.set_xdata(x)
#     plot1.set_ydata(update_y_value)
    
#     figure.canvas.draw()
#     figure.canvas.flush_events()
#     time.sleep(0.1)


# Display

#%%
# https://pythonguides.com/matplotlib-update-plot-in-loop/

# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation

# # Create figure and subplot

# figure, ax = plt.subplots(figsize=(4,3))

# x = []
# y = []

# # Plot

# plot_1, = ax.plot(x, y)
# plt.axis([0, 30, -2, 2])

# # Animation Function

# def animate_plot (i):
#     x = np.linspace(0, 30, 100)
#     y = np.sin(x*i)
#     plot_1.set_data(x, y)
#     return plot_1,

# # Animated Function

# ani = FuncAnimation(figure,
#                     animate_plot,
#                     frames=3,
#                     interval=10)

# # Save as gif

# ani.save('ag_001.gif', fps=10)

# # Display

# plt.show()

  