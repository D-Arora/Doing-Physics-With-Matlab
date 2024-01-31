# animate_001.py
# Ian Cooper
# Jan 2024
# Doing Physics with Python
#   https://d-arora.github.io/Doing-Physics-With-Matlab/
# ANIMATIONS: TRAVELLINFG WAVRS
#    y = asin(kx + wt)

#https://matplotlib.org/stable/api/_as_gen/matplotlib.animation.PillowWriter.html

#https://matplotlib.org/stable/users/explain/animations/animations.html

#https://ajaytech.co/animation/



#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation  
import time
from numpy import pi, sin, cos, linspace


 

#%% INPUTS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# X grid
# Nx = 99; x1= 0; x2 = 50
# x = linspace(x1,x2,Nx)
# # Time span
# Nt = 25; t1 = 0; t2 = 10
# t = linspace(t1,t2,Nt)
# # Wave parameters:
# # amplitude,period, wavelength, frequency, angular frequency, wave number
# A = 5
# T = 2
# L = 10

# f = 1/T; w = 2*pi*f; k = 2*pi/L


#%%
# Initializing a figure in which the graph will be plotted 


fig1 = plt.figure()
ax1 = plt.subplot(111)
line1,  = ax1.plot([],[])
ax1.set_xlim(0,2)
ax1.set_ylim(-1.5,1.5)
x = np.linspace(0,2,1000)
 
def init():
    ax1.set_xlabel('x')
    ax1.set_ylabel('sin(x)')
    ax1.set_title('sine curve')
    line1.set_data([],[])
    return line1,
 
def animate(i):
    y = np.sin(2 * np.pi * (x -  0.01*i))
    line1.set_data(x,y)
   # ax1.set_title('t = %2.3 f ' % i)
    return line1,
 
anim1 = FuncAnimation(fig1,animate,init_func=init,frames=100,interval=100)
 
plt.show()

#fig = plt.figure()  
# font1 = {'family':'Tahoma','color':'black','size':12}
# plt.rcParams['font.family'] = ['Tahoma']
# plt.rcParams['font.size'] = 12

# fig = plt.figure(figsize=(4,4)) 
# ax = fig.add_subplot(111, autoscale_on=False, xlim=(x1,x2), ylim=(-A-1,A+1))
# fig.subplots_adjust(top=0.8, bottom = 0.20, left = 0.18,\
#                     right = 0.93, hspace = 0.20,wspace=0.2)
# #axis = plt.axes(xlim =(x1, x2), ylim =(-A-1, A+1))  
# plt.grid('visible') 
# plt.xlabel('x  [m]',fontsize = 12) 
# plt.ylabel('y [m]',fontsize = 12)
# plt.title('t = %2.3f  s' % t[0])
# #plt.tight_layout()
# # initializing a line variable 
# line, = ax.plot([], [],'b', lw = 3)  
# #time_template = 'time = %.1fs'
# #time_text = ax.text(5,4, '', transform=ax.transAxes)

# def init():  
#     line.set_data([], []) 
#     return line, 
   
# def animate(n):
#      z = np.arange(0,Nt-1,1)
#      y = A*sin(k*x - w*t[n])
#    # ux = [0,np.sin(w*t[0:n])]
#    # vy = [0,np.cos(w*t[0:n])]
#    # u = np.sin(w*t[0:n])
#    # v = np.cos(w*t[0:n])
#     #plt.title('tx = %2.3f  s' % z)
#     # z = np.cos(w*t[n])
#      #x = np.linspace(0, 4, 1000) 
#      u = x
#      v = y
   
#      z = t[n]
#      # plots a sine graph 
#     #y = np.sin(2 * np.pi * (x - 0.01 * i)) 
#      line.set_data(u, v) 
     
#     # ax.set_title('t')
                 
#      #time_text.set_text(time_template % (t[n]))
#      time.sleep(0.1)
      
#      return   line,

#anim = FuncAnimation(fig1, animate, init_func = init, 
#                     frames = 100, interval = 50, blit = True, repeat = False)