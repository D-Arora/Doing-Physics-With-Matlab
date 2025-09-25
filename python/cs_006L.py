# cs_006L.py
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


# LOAD DATA
#[t,theta1] = np.loadtxt('aaa.txt')
#theta2 = np.loadtxt('aaa1.txt')
[t1,theta1] = np.loadtxt('a1.txt')
[t2,theta2] = np.loadtxt('a2.txt')
# Time span
# N = 999         # Time interval t from t1 to t2
# t1 = 0
# t2 = 20
# t = np.linspace(t1,t2,N)

# #%%
# GRAPHICS
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

# Fig 1  t vs theta  ---------------------------------------------------------
fig = plt.figure(figsize = (4, 3))

fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                      right = 0.96, hspace = 0.2,wspace=0.2)

xP = t1; yP = theta1
plt.plot(xP,yP,linewidth=2,color='b')
xP = t2; yP = theta2
plt.plot(xP,yP,linewidth=2,color='r')

# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
plt.grid('visible')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'$\theta$  / $\pi$', fontdict = font1)
plt.savefig('a1.png')

# Fig 2  t vs theta  ---------------------------------------------------------
fig2 = plt.figure(figsize = (4, 3))

fig2.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                      right = 0.96, hspace = 0.2,wspace=0.2)

xP = t1; yP = (abs(theta1-theta2))
plt.plot(xP,yP,linewidth=2,color='b')

# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
plt.grid('visible')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'|$\theta_2$ - $\theta_1$|  / $\pi$', fontdict = font1)
plt.savefig('a2.png')

# Fig 3  t vs theta  ---------------------------------------------------------
fig3 = plt.figure(figsize = (4, 3))

fig3.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                      right = 0.96, hspace = 0.2,wspace=0.2)

xP = t1; yP = np.log(abs(theta1-theta2))
plt.plot(xP,yP,linewidth=2,color='b')

# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
plt.grid('visible')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'log(|$\theta_2$ - $\theta_1$)', fontdict = font1)
plt.savefig('a3.png')

