# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 14:56:59 2024

@author: Owner
"""

#https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.PillowWriter.html

#https://matplotlib.org/stable/users/explain/animations/animations.html



from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation  
import time


N = 99
w = 1
t = np.linspace(0,6,N)
dt = t[2]-t[1]
x1 = np.sin(w*t)
y = np.zeros(N)
z = np.zeros(N)
y1 = np.cos(w*t)
# initializing a figure in which the graph will be plotted 
fig = plt.figure()  
   
# marking the x-axis and y-axis 
axis = plt.axes(xlim =(-1, 1),  
                ylim =(-1, 1))  
  
# initializing a line variable 
line, = axis.plot([], [],'o', lw = 3)  


def init():  
    line.set_data([], []) 
    return line, 
   
def animate(n):
   # t = np.linspace(0,10,N)
    z = np.arange(0,N-1,1)
    #  x = t
   # ux = [0,np.sin(w*t[0:n])]
   # vy = [0,np.cos(w*t[0:n])]
   # u = np.sin(w*t[0:n])
   # v = np.cos(w*t[0:n])
    
    # z = np.cos(w*t[n])
     #x = np.linspace(0, 4, 1000) 
    u = x1[0:n]
    v = y1[0:n]
   
     # plots a sine graph 
    #y = np.sin(2 * np.pi * (x - 0.01 * i)) 
    line.set_data(u, v) 
    time.sleep(0.1)
      
    return line,  

anim = FuncAnimation(fig, animate, init_func = init, 
                     frames = 100, interval = 50, blit = True, repeat = False)