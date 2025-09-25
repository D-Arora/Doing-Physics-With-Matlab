# cs_006.py
# Ian Cooper
# Jan 2024
# COMPLEX SYSTEMS
#  TIME DEPENDENT DYNAMICAL SYSTEMS
#      	PENDULUM: free, damped, forced motion 


# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_006.pdf

# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
import math
from numpy import pi, sin, cos
from matplotlib.animation import FuncAnimation, PillowWriter  

tStart = time.time()

#%%
# FREE, DAMPED AND DRIVEN MOTION OF A SIMPLE PENDULUM
#  SI units
#  angular displacement theta  [rad]
#  angular velocity (angular frequency) omega  [rad/s]

# INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Initial angular displacment  [ -360 < theta(0) < +360  deg]
thetaDeg0 = 90

# initial angular velocity  [rad/s]
omega0 = 0

# Acceleration due to gravity
g = 9.8     
# length of pendulum  [m]   --> period = 1s
L = 3*g/(4*pi**2)

# Damping coefficient  0 <= b < 10
b = 0   
        
# Driving force
# Amplitude of driving force    0 <= AD <= 10
AD = 10
# frequency of driving force   fD ~ 1 Hz                   
fD = 1.2

# Time span
N = 999         # Time interval t from t1 to t2
t1 = 0
t2 = 50


#%%
# SETUP:Free pendulum:  Period  [s] / Frequency  /  Angular frequency
T0 = 2*np.pi*(L/g)**0.5 
f0 = 1/T0  
w0 = 2*pi/T0

# Angular frequency of driving force and period
wD = 2*pi*fD
TD = 1/fD

  
#%%  SOLVE ODE
# Solve ODE for x,y    x = theta   y = omega 
def lorenz(t, state):    
    x, y = state
    dx = y
    dy = -(g/L)*sin(x) - b*y + AD*sin(wD*t)
    return [dx, dy]  

tSpan = np.linspace(t1,t2,N)
t = tSpan
theta0 = thetaDeg0*pi/180            # degrees to radians
u0 = [theta0, omega0]
sol = odeint(lorenz, u0, tSpan, tfirst=True)
theta = sol[:,0]/pi     # angular displacement [rad/pi]
omega = sol[:,1]        # angular velocity  [rad/s]      


#%%

# Wrapping theta
theta = np.arctan2( np.sin(theta), np.cos(theta))



      
#%%
# GRAPHICS
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

# Fig 1  theta (x) vs t  ---------------------------------------------------------
fig1 = plt.figure(figsize = (4, 3))

fig1.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                      right = 0.96, hspace = 0.2,wspace=0.2)

xP = t; yP = theta
plt.plot(xP,yP,linewidth=2,color='b')
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
plt.grid('visible')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'$\theta$  / $\pi$', fontdict = font1)
# plt.savefig('cs001.png')



#%%
x1 = -1; x2 = 1; y1= -1; y2 = 1

#  # Initializing a figure in which the graph will be plotted 
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize=(4,4)) 
ax = fig.add_subplot(111, autoscale_on=False, xlim=(x1,x2), ylim=(y1,y2))
fig.subplots_adjust(top=0.98, bottom = 0.20, left = 0.10,\
                    right = 0.90, hspace = 0.20,wspace=0.2)

ax.set_aspect('equal', adjustable='box') 
ax.set_xticks([])
ax.set_yticks([])   
#plt.grid('visible') 
plt.xlabel('x',fontsize = 12) 
plt.ylabel('y',fontsize = 12)
time_text = ax.text(-0.8,0.8, '')
# initializing a line variable 
line, = ax.plot([], [],'bo',markersize = 10) 
line1, = ax.plot([], [],'ro', markersize = 8)  
line2, = ax.plot([], [],'r-',lw = 3)
#%% FUNCTIONS
def init():  
    line.set_data([], [])
    time_text.set_text('')
    
    return line, 
   
def animate(n):
     # u = L*sin(omega[0:n]*t[0:n])
     # v = -L*cos(omega[0:n]*t[0:n])
      u = L*sin(omega[n]*t[n])
      v = -L*cos(omega[n]*t[n])
        
      line.set_data([u], [v]) 
      line1.set_data([0], [0])
      line2.set_data([0,u], [0,v])
      
      time_text.set_text('time = %.2f' % t[n] + ' s')  
      time.sleep(0.1)
      return   line, line1, line2, time_text,

anim = FuncAnimation(fig, animate, init_func = init, 
                      frames = 199, interval = 20, blit = True, repeat = False)

# # anim.save('ag_001.gif', fps = 2)   



