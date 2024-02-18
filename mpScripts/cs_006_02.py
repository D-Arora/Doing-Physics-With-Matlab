# cs_006_01.py
# Ian Cooper
# Jan 2024
# COMPLEX SYSTEMS
#  TIME DEPENDENT DYNAMICAL SYSTEMS
#      	PENDULUM: free, damped, forced motion 
#       CHAOTIC DYNAMICS

# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_006.pdf


# https://scipython.com/book/chapter-6-numpy/examples/finding-a-best-fit-straight-line/


# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
import math
from numpy import pi, sin, cos
Polynomial = np.polynomial.Polynomial


tStart = time.time()

#%%
# FREE, DAMPED AND DRIVEN MOTION OF A SIMPLE PENDULUM
#  SI units
#  angular displacement theta  [rad]
#  angular velocity (angular frequency) omega  [rad/s]

# INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Initial angular displacment [radian]  
theta10 = 0
theta20 = 0.001    
# initial angular velocity  [rad/s]
omega0 = 0

# Driving force: Amplitude / drive strength (gamma)  / Frequency
AD = 1.503

fD = 1.0   

# Time span: t1 to t2  
N = 9999          
t1 = 0.0

t2 = 40


#%%  SETUP
# DRIVE FORCE
wD = 2*pi*fD
TD = 1/fD
# Natural frequency and period
w0 = 1.5*wD
f0 = w0/(2*pi) 
T0 = 1/f0

# Damping constant  >>>>
b = w0/4

# Pendulum
g = 9.8
L = g/w0**2

 
#%%  SOLVE ODE
# Solve ODE for x,y    x = theta   y = omega 
def lorenz(t, state):    
    x, y = state
    dx = y
    dy = -w0**2*sin(x) - 2*b*y + AD*w0**2*cos(wD*t)
    return [dx, dy]  

tSpan = np.linspace(t1,t2,N)
t = tSpan
dt = t[1] - t[0]
u0 = [theta10, omega0]
sol = odeint(lorenz, u0, tSpan, tfirst=True)
theta1 = sol[:,0]     # angular displacement [rad/pi]
omega1 = sol[:,1]        # angular velocity  [rad/s]      

u0 = [theta20, omega0]
sol = odeint(lorenz, u0, tSpan, tfirst=True)
theta2 = sol[:,0]     # angular displacement [rad/pi]
omega2 = sol[:,1]        # angular velocity  [rad/s]    

#theta1 = theta1%(2*pi)
#theta2 = theta2%(2*pi)
delta = abs(theta2 - theta1)

#%%
# GRAPHICS
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12


#%% Fig 0   Time evolution theta2 and theta1
plt.rcParams["figure.figsize"] = (6,2.7)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.89, bottom = 0.21, left = 0.150,\
                    right = 0.96, hspace = 0.36,wspace=0.50)

axes.set_ylabel(r'$\theta$ / $\pi$',color= 'black',fontsize = 12)
axes.set_xlabel('$t$  [ s ]',color = 'black',fontsize = 12)
axes.set_title('$\gamma = $ %2.4f ' % AD )
#axes.set_xlim([20, 40])
#axes[R].set_ylim([0,1])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_yticks(np.arange(0,1.1,0.2))
axes.xaxis.grid()
axes.yaxis.grid()

xP = t; yP = theta1/pi
axes.plot(xP, yP, 'b',lw = 3)
yP = theta2/pi
axes.plot(xP, yP, 'r',lw = 1)

plt.savefig('a0.png') 


#%% Fig 2   Time evolution delta
plt.rcParams["figure.figsize"] = (6,2.7)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.89, bottom = 0.21, left = 0.150,\
                    right = 0.96, hspace = 0.36,wspace=0.50)

axes.set_ylabel(r'| $\Delta$ $\theta$ / $\pi$ | ',color= 'black',fontsize = 12)
axes.set_xlabel('$t$  [ s ]',color = 'black',fontsize = 12)
axes.set_title('$\gamma = $ %2.4f ' % AD )
     
#axes.set_xlim([20, 40])
#axes[R].set_ylim([0,1])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_yticks(np.arange(0,1.1,0.2))
axes.xaxis.grid()
axes.yaxis.grid()

xP = t; yP = delta  
axes.plot(xP, yP, 'b')

plt.savefig('a1.png') 


#%% Fig 2   Time evolution delta
plt.rcParams["figure.figsize"] = (6,2.7)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.89, bottom = 0.21, left = 0.190,\
                    right = 0.96, hspace = 0.36,wspace=0.50)

axes.set_ylabel(r'| $\Delta$ $\theta$ / $\pi$ | ',color= 'black',fontsize = 12)
axes.set_xlabel('$t$  [ s ]',color = 'black',fontsize = 12)
axes.set_title('$\gamma = $ %2.4f ' % AD )
#axes.set_xlim([20, 40])
#axes[R].set_ylim([0,1])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_yticks(np.arange(0,1.1,0.2))
axes.xaxis.grid()
axes.yaxis.grid()
plt.yscale("log")  
 
xP = t; yP = delta  
axes.plot(xP, yP, 'b')

plt.savefig('a2.png') 


#%%
tEnd = time.time()
tE = tEnd - tStart
print('\nExecution time =  %2.3f s' % tE )


