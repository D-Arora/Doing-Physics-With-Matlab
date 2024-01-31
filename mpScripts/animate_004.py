# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 15:06:19 2024

@author: Owner
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:41:09 2020

Simulate then animate the simple pendulum motion.
Method:  Runge-Kutta integration with Scipy

@author: Prof. Roy, University of Delaware
"""

# https://research.me.udel.edu/~vroy/python/tuto7/


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import scipy.integrate as integrate
########## parameters ####################
g = 9.81 # [m/s^2]
l = 0.5 # [m]
gamma = 0. # [1/s] damping coefficient
omega0 = np.sqrt(g/l)

############# initial conditions  ################

theta= 170 # [degree]  # initial angle
theta = theta*np.pi/180
omega = 0  # initial angular velocity

# given the initial condition, we can compute the period of oscillations
# by a quadrature. use the result to set up max time of simulation

def f(x):
    return 1/np.sqrt(np.cos(x)-np.cos(theta))

Tper = integrate.quad(f ,0,theta )
T0 = 2.0*np.sqrt(2.0)*Tper[0]/omega0
Tmax = 2*T0 # set Tmax as a multiple of Tper
print('Tmax=',Tmax)

# Graphical and simulation parameters
 
xmax =  1.05*l
ymax =  1.05*l
Nstep= 800   # number of time steps /number of frames
             # increase value of Nstep for longer simulation
dt= Tmax/Nstep  # time increment
d=0.05 # size of fixed support

# Set up the figure, the axis, and the line object we want to animate.
# This object is a line segment connecting two points: the origin (0,0)
# and the variable coordinates of the pendulum end point

fig = plt.figure()
ax = plt.axes(xlim=(-xmax, xmax), ylim=(-ymax, ymax))

line, = ax.plot([], [], 'o', ls='-', ms=8, markevery=[0,-1])
xR,yR = -d,-d # lower left corner of rectangle
dx,dy= 2*d, 2*d  # dimensions of rectangle
rect = plt.Rectangle((xR,yR), dx,dy, facecolor= 'black', edgecolor= 'black')
ax.add_patch(rect)
plt.axis('off')
ax.set_aspect('equal')

# define the equations of motion
def F(X,t):
    [theta, omega] = X
    dthetadt = omega
    domegadt = -g/l* np.sin(theta) - gamma *omega
    return dthetadt, domegadt

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function: compute pendulum position in Cartesian coordinate 
# at each time step by integration of the ODE. We use scipy ODE integrator.
    
def animate(i):
    
    global theta, omega
    timeIn= i*dt
    timeOut= (i+1)*dt
    Xinit = [theta, omega]
    odeout = integrate.odeint(F, Xinit, [timeIn, timeOut])[1]
    theta = odeout[0] 
    omega = odeout[1]
    theta = (theta + np.pi) % (2*np.pi) - np.pi # bring theta back to [-pi, pi]
    x,y = l*np.sin(theta), -l*np.cos(theta)
    line.set_data([0, x], [0,y]) # modify the pendulum's position    
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=Nstep, interval=20, blit=True)

# save the animation as an mp4.  
#anim.save('pendulum_animation_odeint.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
anim.save('ag_010.gif')

plt.show()