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
from numpy import pi, sin

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
omega0 = 12.566
# Acceleration due to gravity
g = 9.8     
# length of pendulum  [m]   --> period = 1s
L = g/(4*pi**2)

# Damping coefficient  0 <= b < 10
b = 0   
        
# Driving force
# Amplitude of driving force    0 <= AD <= 10
AD = 0
# frequency of driving force   fD ~ 1 Hz                   
fD = 1.2566

# Time span
N = 999         # Time interval t from t1 to t2
t1 = 0
t2 = 5


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


#%%   Find peaks in x vs t plot
peaks, _ = find_peaks(theta)
numT = len(peaks)-1
Tpeaks = (t[peaks[-1]] - t[peaks[0]])/numT
fpeaks = 1/Tpeaks

#%% FOURIER TRANSFORM - frequency spectrum

def simpson1d(f,xMin,xMax):
    N = len(f)
    h = (xMax - xMin) / (N - 1)
    
    integral = (h/3) * (f[0] + 2*sum(f[:N-2:2]) \
            + 4*sum(f[1:N-1:2]) + f[N-1])
 
    if N%2 == 0:
        integral = 'N must be an odd number'
        print('integral')
    return integral

p  = 1j*2*pi
Fmin = -2
Fmax = 2
nF = 299
F = np.linspace(Fmin,Fmax,nF)
HR = np.zeros(nF); HI = HR; H = HR; H = HR+HR*1j
for c in range(nF):
     g = theta*np.exp(p*F[c]*t)
     H[c] = simpson1d(g,t1,t2)

psd = 2*np.conj(H)*H   
psd = np.real(psd/max(psd))   

Find = np.where(np.real(psd)>= 0.98)
II = Find[-1]
II = II[-1]
Fpeak = abs(F[II])


#%% sinusoidal motion   SHM
A0 = max(abs(theta))
phi = math.asin(u0[0]/(A0*pi))  
xS = A0*sin(2*pi*f0*t + phi)

  
#%% CONSOLE OUTPUT
print(' ')
print('Input parameters')
print('theta(0) = %2.3f deg' % theta0, '    omega(0) = %2.3f  rad/s' % omega0)
print('b = %2.3f' % b, ' fD = %2.3f  Hz' % fD \
      ,'  AD = %2.3f ' % AD)
print('Results')
print('Free vibration: T0 = %2.3f s' % T0, '   f0 = %2.3f  Hz' % f0)
print('Peaks theta vs t:  Tpeak = %2.3f s'  % Tpeaks, '  fPeaks = %2.3f  Hz' %fpeaks)
print('psd:    fPeak = %2.3f  Hz' % Fpeak)

      
#%%
# GRAPHICS
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

# Fig 1  theta (x) vs t  ---------------------------------------------------------
fig = plt.figure(figsize = (4, 3))

fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                     right = 0.96, hspace = 0.2,wspace=0.2)

xP = t; yP = theta
plt.plot(xP,yP,linewidth=2,color='b')
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
plt.grid('visible')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'$\theta$  / $\pi$', fontdict = font1)
#plt.savefig('cs001.png')

# Fig 2  omega (y) vs t  ---------------------------------------------------------
fig = plt.figure(figsize = (4, 3))

fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                     right = 0.96, hspace = 0.2,wspace=0.2)

xP = t; yP = omega
plt.plot(xP,yP,linewidth=2,color='b')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'$\omega$   [rad/s]', fontdict = font1)

plt.grid('visible')
#plt.savefig('cs002.png')

# Fig 3  Phase portrait  omega vs theta --------------------------------
fig = plt.figure(figsize = (4, 3))

fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                     right = 0.96, hspace = 0.2,wspace=0.2)

xP = theta; yP = omega
plt.plot(xP,yP,linewidth=2,color='b')
plt.ylabel(r'$\omega$   [ rad/s ]', fontdict = font1)
plt.xlabel(r'$\theta$ $\pi$', fontdict = font1)
plt.grid('visible')
# #plt.savefig('cs003.png')

# Fig 4  Frequency spectrum --------------------------------
fig = plt.figure(figsize = (4, 3))

fig.subplots_adjust(top = 0.94, bottom = 0.24, left = 0.22,\
                     right = 0.96, hspace = 0.2,wspace=0.2)

xP = F[F>0]; yP = psd[F>0]
plt.plot(xP,yP,linewidth=2,color='b')
plt.xlabel(r'$\ f   $  [ Hz ]', fontdict = font1)
plt.ylabel(r'$\ psd $', fontdict = font1)
plt.grid('visible')
#plt.savefig('cs004.png')             

#%% SUBPLOTS   Figure 5

fig = plt.figure(figsize = (8.5, 6))
fig.subplots_adjust(top = 0.97, bottom = 0.12, left = 0.10,\
                     right = 0.98, hspace = 0.32,wspace=0.24)

plt.subplot(2, 2, 1)
xP = t; yP = theta
plt.plot(xP,yP,linewidth=2,color='b')
xP = t; yP = xS
plt.plot(xP,yP,linewidth=1,color='r')
plt.grid('visible')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'$\theta$ / $\pi$ ', fontdict = font1)

plt.subplot(2, 2, 2)
xP = t; yP = omega
plt.plot(xP,yP,linewidth=2,color='b')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'$\omega$  [rad/s]', fontdict = font1)
plt.grid('visible')

plt.subplot(2, 2, 3)
xP = theta; yP = omega
plt.plot(xP,yP,linewidth=2,color='b')
plt.ylabel(r'$\omega$   [ rad/s ]', fontdict = font1)
plt.xlabel(r'$\theta$ / $\pi$', fontdict = font1)
plt.grid('visible')

plt.subplot(2, 2, 4)
xP = F[F>0]; yP = psd[F>0]
plt.plot(xP,yP,linewidth=2,color='b')
plt.xlabel(r'$\ f   $  [ Hz ]', fontdict = font1)
plt.ylabel(r'$\ psd $', fontdict = font1)
plt.grid('visible')

plt.savefig('cs005.png')  

tEnd = time.time()
tRun = tEnd - tStart
print('  ')
print(tRun)




