# cs_006_01.py
# Ian Cooper
# Jan 2024
# COMPLEX SYSTEMS
#  TIME DEPENDENT DYNAMICAL SYSTEMS
#      	PENDULUM: free, damped, forced motion 
#       CHAOTIC DYNAMICS

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

tStart = time.time()

#%%
# FREE, DAMPED AND DRIVEN MOTION OF A SIMPLE PENDULUM
#  SI units
#  angular displacement theta  [rad]
#  angular velocity (angular frequency) omega  [rad/s]

# INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Initial angular displacment [radian]  
theta0 = -pi/2       
# initial angular velocity  [rad/s]
omega0 = 0

# Driving force: Amplitude / drive strength (gamma)  / Frequency
AD = 1.06

fD = 1.0   


NS = -1000

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


#%%
# GRAPHICS
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize = (5, 4))

fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.25,\
                      right = 0.96, hspace = 0.2,wspace=0.2)
plt.grid('visible')
plt.xlabel(r'$\gamma$  ', fontdict = font1)
plt.ylabel(r'$\omega$  [ rad/s ]', fontdict = font1)    
 
#%%  SOLVE ODE
# Solve ODE for x,y    x = theta   y = omega 
def lorenz(t, state):    
    x, y = state
    dx = y
    dy = -w0**2*sin(x) - 2*b*y + AD*w0**2*cos(wD*t)
    return [dx, dy]  

#%%
N = 1000
dt = 1/N
tMax = 600
tMin = 0
t = np.arange(tMin,tMax,dt)

ind = N*np.arange(500,600,1)
ind = ind.astype(int)
tP = t[ind]
u0 = [theta0, omega0]
Ng = 600
#gamma = np.linspace(1.061,1.087,Ng)
gamma = np.linspace(1.03,1.53,Ng)


#%%
for c in range(Ng-1):
    AD = gamma[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    theta = sol[:,0]     # angular displacement [rad/pi]
    omega = sol[:,1]        # angular velocity  [rad/s]      
    
    xP = np.ones(len(ind))*gamma[c]
    yP = omega[ind]
    plt.plot(xP,yP,'bo', ms = 0.5)




#%%
tEnd = time.time()
tE = (tEnd - tStart)/60
print('\nExecution time =  %2.3f min' % tE )

# #%%

# N = 10
# dt = 1/N
# tMax = 60
# tMin = 0

# t = np.arange(tMin,tMax,dt)

# plt.savefig('a0.png') 