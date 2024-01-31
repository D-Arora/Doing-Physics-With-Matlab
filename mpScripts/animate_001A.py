# animate_001.py
# Ian Cooper
# Jan 2024
# Doing Physics with Python
#   https://d-arora.github.io/Doing-Physics-With-Matlab/

# ANIMATIONS: TRAVELLING WAVES
#    y = Asin(kx +/- wt)

#https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.PillowWriter.html

#https://matplotlib.org/stable/users/explain/animations/animations.html

#https://research.me.udel.edu/~vroy/python/tuto5/


#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace


#%% INPUTS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# X grid
Nx = 299; x1= 0; x2 = 50
x = linspace(x1,x2,Nx)
# Time span
Nt = 125; t1 = 0; t2 = 10
t = linspace(t1,t2,Nt)
# Wave parameters:
# amplitude,period, wavelength, frequency, angular frequency, wave number
A = 5
T = 2
L = 20

f = 1/T; w = 2*pi*f; k = 2*pi/L


#%%  # Initializing a figure in which the graph will be plotted 
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize=(4,4)) 
ax = fig.add_subplot(111, autoscale_on=False, xlim=(x1,x2), ylim=(-A-1,A+2))
fig.subplots_adjust(top=0.98, bottom = 0.20, left = 0.18,\
                    right = 0.93, hspace = 0.20,wspace=0.2)
plt.grid('visible') 
plt.xlabel('x  [m]',fontsize = 12) 
plt.ylabel('y [m]',fontsize = 12)
time_text = ax.text(14,6.2, '')
# initializing a line variable 
line,  = ax.plot([], [],'b', lw = 3)  
line1, = ax.plot([], [],'ro', ms = 8)  

#%% FUNCTIONS
def init():  
    line.set_data([], [])
    line1.set_data([], [])
   # time_text.set_text('')
    
    return line, line1 
   
def animate(n):
     y = A*sin(k*x - w*t[n])
     
     y1 = A*sin(k*x[50] - w*t[n])
     
     u = x
     v = y
     u1 = x[50]
     v1 = y1
   
     line.set_data([u], [v]) 
     line1.set_data([u1], [v1]) 
     
     time_text.set_text('time = %.2f' % n + ' s')  
     time.sleep(0.1)
     return   line, line1, time_text,

anim = FuncAnimation(fig, animate, init_func = init, 
                     frames = 101, interval = 20, blit = True, repeat = False)

# anim.save('ag_001.gif', fps = 2)   



