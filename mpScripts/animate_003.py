# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 14:21:20 2024

@author: Owner
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:41:09 2020

Animation of projectile motion with time clock.

@author: Prof Roy
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation


########## parameters ####################
v0 = 50 # [m/s]
g = 9.81 # [m/s^2]
alpha= 50 # [deg]
D= 220 # [m]
d= 10 # [m]
alpha = alpha*np.pi/180 # launch angle [rad]
Tmax = D/(v0*np.cos(alpha)) 
xmax =  D+2*d
ymax =  0.55*(v0*np.sin(alpha))**2 /g 
Nstep= 400   # number of time steps /number of frames
dt= Tmax/Nstep  # time increment
Vscale=0.2
h= D*np.tan(alpha)-0.5*g*(D/v0/np.cos(alpha))**2

# First set up the figure, the axis, and the line object we want to animate

fig = plt.figure()
ax = plt.axes(xlim=(0, xmax), ylim=(0, ymax))
plt.xlabel('x(m): horizontal distance')
plt.ylabel('y(m): height')
time_text = ax.text(0.1*xmax, 0.9*ymax, '')

line, = ax.plot([], [], 'o', ls='-', ms=8, markevery=[0,-1])
xR,yR = D-d,0. # lower left corner of rectangle
dx,dy= 2*d, h  # dimensions of rectangle
rect = plt.Rectangle((xR,yR), dx,dy, facecolor= 'black', edgecolor= 'black')
ax.add_patch(rect)


# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line,

# animation function: compute projectile position in Cartesian coordinate vs time
def animate(i):
    ttmax=i*dt
    t= np.linspace(0, ttmax, 2*i)
    x = v0*np.cos(alpha)*t
    y = v0*np.sin(alpha)*t -0.5*g*t**2
    line.set_data(x, y)
    time_text.set_text('time = %.2f' % ttmax + ' s')
    return line, time_text,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Nstep, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  
anim.save('projectile_animation_clocked.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
plt.show()