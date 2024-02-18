# cs_006.py
# Ian Cooper
# Jan 2024
# COMPLEX SYSTEMS
#  TIME DEPENDENT DYNAMICAL SYSTEMS
#      	PENDULUM: free, damped, forced motion 


# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_006.pdf


# https://www.evernote.com/shard/s140/client/snv?isnewsnv=true&noteGuid=0724815b-79a9-4357-9e85-416c33cb1b69&noteKey=e2b0667446e6f7d74181969ed0c7c357&sn=https%3A%2F%2Fwww.evernote.com%2Fshard%2Fs140%2Fsh%2F0724815b-79a9-4357-9e85-416c33cb1b69%2Fe2b0667446e6f7d74181969ed0c7c357&title=Chapter%2B3%2BOscillatory%2BMotion%2Band%2BChaos

# https://www.evernote.com/shard/s140/client/snv?isnewsnv=true&noteGuid=0724815b-79a9-4357-9e85-416c33cb1b69&noteKey=e2b0667446e6f7d74181969ed0c7c357&sn=https%3A%2F%2Fwww.evernote.com%2Fshard%2Fs140%2Fsh%2F0724815b-79a9-4357-9e85-416c33cb1b69%2Fe2b0667446e6f7d74181969ed0c7c357&title=Chapter%2B3%2BOscillatory%2BMotion%2Band%2BChaos

#
# https://physics.stackexchange.com/questions/398995/issue-with-bifurcation-plot-for-driven-pendulum


# *** 
#https://www.physics.purdue.edu/~hisao/book/www/Computational%20Physics%20using%20MATLAB%20v2.pdf



# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
import math
from numpy import pi, sin, cos
import array

tStart = time.time()

#%%
# FREE, DAMPED AND DRIVEN MOTION OF A SIMPLE PENDULUM
#  SI units
#  angular displacement theta  [rad]
#  angular velocity (angular frequency) omega  [rad/s]

# INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Initial angular displacment  [ -360 < theta(0) < +360  deg]
thetaDeg0 = 0

# initial angular velocity  [rad/s]
omega0 = 0
# Acceleration due to gravity
g = 9.8     
# length of pendulum  [m]   --> period = 1s
L = g/(4*pi**2)
L = g
# Damping coefficient  0 <= b < 10
b = 0.5

        
# Driving force
# Amplitude of driving force    0 <= AD <= 10
AD = 1.35
# frequency of driving force   fD ~ 1 Hz                   
fD = (2/3)/(2*pi) #0.71

# Time span
N = 5000         # Time interval t from t1 to t2
t1 = 0.0
t2 = 100


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
    # dy = -(g/L)*sin(x) - b*y + AD*cos(wD*t)
    dy = -9*pi**2*sin(x) - (1.5*pi)*y + AD*cos(2*pi*t)
    return [dx, dy]  

tSpan = np.linspace(t1,t2,N)
t = tSpan
theta0 = thetaDeg0*pi/180            # degrees to radians
u0 = [theta0, omega0]


#  Driving frequency span  f
#f1 = 0.5; f2 = 5; nF = 99
#f = np.linspace(f1,f2,nF)
#w = 2*pi*f
nS = 10000; A1 = 0.6; A2 = 2.0
A = np.linspace(A1,A2,nS+1)

# GRAPHICS
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize = (6, 4))

fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.15,\
                      right = 0.96, hspace = 0.2,wspace=0.2)

plt.xlim([A1,A2])
plt.grid('visible') 
plt.ylim([-1.1,1.1])

plt.xlabel(r'strength of drivng force  $A_D$ ', fontdict = font1)
plt.ylabel(r'$\theta /  \pi$', fontdict = font1)

theta = np.zeros(nS+1)

for c in range(nS):  #nF
    AD = A[c]*9*pi**2
    sol = odeint(lorenz, u0, tSpan, tfirst=True)
    #theta = sol[:,0]/pi     # angular displacement [rad/pi]
    omega = sol[:,1]        # angular velocity  [rad/s]      
   # omegaR = omega[omega>0]
   # omegaR = omegaR[-100:-10]
   # y = omegaR
    y = sol[:,0]
# Findpeaks
    
   
    theta[c] = y[-1]
    
    while theta[c] > +pi:
         theta[c] = theta[c] - pi
    while theta[c] < -pi:
         theta[c] = theta[c] + pi
   # ind = indexMax[indexMax > 0]
   # ind = ind.astype(int)

    # OR1 = max(omegaR); OR2 = OR1/2
    # #ind = find_peaks(omegaR,distance =  20, height = (OR2,OR1), threshold = OR2) [0]
    # ind = find_peaks(omegaR)[0]
   # Lind = len(ind)
   # wPeaks = np.zeros(Lind)
   # wPeaks = y[ind]
   # fPeaks = np.ones(Lind)*wD/(2*pi)
   # xP = fPeaks; yP = wPeaks
   # plt.plot(xP,yP,'bo',markersize = 1)



xP = A; yP = theta/pi
plt.plot(xP,yP,'bo',markersize = 1)    
    
#%%
# Wrapping theta
#theta = np.arctan2( np.sin(theta), np.cos(theta))




# #%% CONSOLE OUTPUT
# print(' ')
# print('Input parameters')
# print('theta(0) = %2.3f deg' % theta0, '    omega(0) = %2.4f  rad/s' % omega0)
# print('b = %2.3f' % b, ' fD = %2.3f  Hz' % fD \
#       ,'  AD = %2.3f ' % AD)
# print('Results')
# print('Free vibration: T0 = %2.3f s' % T0, '   f0 = %2.3f  Hz' % f0)


      
# #%%
# # GRAPHICS
# font1 = {'family':'Tahoma','color':'black','size':12}
# plt.rcParams['font.family'] = ['Tahoma']
# plt.rcParams['font.size'] = 12



# #%% Fig 2  t vs omega  ---------------------------------------------------------
# fig = plt.figure(figsize = (4, 3))

# fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
#                      right = 0.96, hspace = 0.2,wspace=0.2)

# xP = t; yP = omega
# plt.plot(xP,yP,linewidth=1,color='b')
# plt.xlabel(r'$t$   [ s ]', fontdict = font1)
# plt.ylabel(r'$\omega$   [rad/s]', fontdict = font1)

# plt.grid('visible')
# #plt.savefig('cs002.png')



#%%
tEnd = time.time()
tRun = tEnd - tStart
print('  ')
print('Execution time = %2.2f  s' % tRun)

#%%
# fig, ax = plt.subplots(1)
# ind = find_peaks(omega>0,distance =  20)[0]
# Lind = len(ind)
# z = np.arange(0,Lind,1)
# plt.plot(z,omega[ind],'o')



