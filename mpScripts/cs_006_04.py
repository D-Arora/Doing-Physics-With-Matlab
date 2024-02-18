# cs_006_01.py
# Ian Cooper           Feb 2024
# COMPLEX SYSTEMS
#  TIME DEPENDENT DYNAMICAL SYSTEMS
#      	PENDULUM: free, damped, forced motion 
#       CHAOTIC DYNAMICS

# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_006A.htm


# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos

tStart = time.time()

#%%   SETUP
#     FREE, DAMPED AND DRIVEN MOTION OF A SIMPLE PENDULUM
#     SI units
#     angular displacement theta  [rad]
#     angular velocity (angular frequency) omega  [rad/s]

# Initial angular displacment [radian]  
theta0 = -pi/2      
# initial angular velocity  [rad/s]
omega0 = 0
# Driving force: Amplitude / drive strength (gamma)  / Frequency
gamma = 1.073
fD = 1.0   
wD = 2*pi*fD
TD = 1/fD
# Natural frequency and period
w0 = 1.5*wD
f0 = w0/(2*pi) 
T0 = 1/f0

# Damping constant  >>>>
b = w0/4  # 4


 
#%%  SOLVE ODE
# Solve ODE for x,y    x = theta   y = omega 
def lorenz(t, state):    
    x, y = state
    dx = y
    dy = -w0**2*sin(x) - 2*b*y + gamma*w0**2*cos(wD*t)
    return [dx, dy]  


N = 1000
dt = 1/N
tMax = 500
tMin = 0


t = np.arange(tMin,tMax,dt)
t1 = 400
t2 = tMax
ind = N*np.arange(t1,t2,1)
ind = ind.astype(int)
tP = t[ind]


u0 = [theta0, omega0]
sol = odeint(lorenz, u0, t, tfirst=True)
theta = sol[:,0]     # angular displacement [rad/pi]
omega = sol[:,1]        # angular velocity  [rad/s]      
 

#%%
# Wrapping theta
#thetaW = np.arctan2( np.sin(theta), np.cos(theta))
# theta  -pi to +pi

NP = len(t)
thetaP = np.zeros(NP)
for c in range(NP):
     thetaP[c] = theta[c]
     while thetaP[c] > pi:
           thetaP[c] = thetaP[c] - 2*pi
     while thetaP[c] < -pi:
          thetaP[c] = thetaP[c] + 2*pi
        
        


  
#%% CONSOLE OUTPUT
print(' ')
p = theta[0]/pi
print('Initial conditions')
print('  theta(0)/pi = %2.3f ' % p, '  omega(0) = %2.4f rad/s' % omega0 \
     ,'  Damping: b = %2.3f' % b )
print('Free vibration')
print('               T0 = %2.3f s' % T0, '  f0 = %2.3f Hz' % f0,  '   w0 = %2.3f rad/s' % w0 )

print('Driving force')
print('  AD = %2.3f' % gamma, '  TD = %2.3f s' % TD \
      ,'  fD = %2.3f Hz' % fD, '   wD = %2.3f rad/s' %wD)


      
#%% GRAPHICS   POINCARE SECTION
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

plt.rcParams["figure.figsize"] = (4,3)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.91, bottom = 0.18, left = 0.180,\
                    right = 0.96, hspace = 0.36,wspace=0.50)


axes.set_ylabel(r'$\omega$ [ rad/s ]',color= 'black',fontsize = 12)
axes.set_xlabel(r'$\theta$ / $\pi$ ',color = 'black',fontsize = 12)
axes.set_title(r'$\gamma = $ %2.4f ' % gamma )
     
axes.set_xlim([-0.15, -0.02])
axes.set_ylim([14, 20])
#axes.set_ylim([1.2*min(omega[ind]),1.2*max(omega[ind])])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_yticks(np.arange(0,1.1,0.2))
axes.xaxis.grid()
axes.yaxis.grid()
xP = thetaP[ind]/pi; yP = omega[ind]
axes.plot(xP, yP, 'bo',ms = 8)

plt.savefig('a1.png') 


plt.rcParams["figure.figsize"] = (4,3)
fig2, axes = plt.subplots(nrows=1, ncols=1)
fig2.subplots_adjust(top = 0.91, bottom = 0.18, left = 0.180,\
                    right = 0.96, hspace = 0.36,wspace=0.50)


axes.set_ylabel(r'$\omega$ [ rad / s ]',color= 'black',fontsize = 12)
axes.set_xlabel(r'$\theta$ / $\pi$ ',color = 'black',fontsize = 12)
axes.set_title(r'$\gamma = $ %2.4f ' % gamma )
     
#axes.set_xlim([-1, 1])
#axes.set_ylim([0, 2])
#axes.set_ylim([1.2*min(omega[ind]),1.2*max(omega[ind])])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_yticks(np.arange(0,1.1,0.2))
axes.xaxis.grid()
axes.yaxis.grid()
xP = theta/pi; yP = omega
axes.plot(xP, yP, 'b',lw = 1)
xP = theta[-10000:-1]/pi; yP = omega[-10000:-1]
axes.plot(xP, yP, 'r',lw = 2)
xP = theta[0]/pi; yP = omega[0]
axes.plot(xP, yP, 'og',ms = 8)

plt.savefig('a2.png') 

#%%
tEnd = time.time()
tE = (tEnd - tStart)/60


print('\nExecution time =  %2.3f s' % tE )


