# cs_006_01.py
# Ian Cooper     Feb 2024
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
# FREE, DAMPED AND DRIVEN MOTION OF A SIMPLE PENDULUM
#  SI units
#  angular displacement theta  [rad]
#  angular velocity (angular frequency) omega  [rad/s]
# Initial angular displacment [radian]  
theta0 = 0.     
# initial angular velocity  [rad/s]
omega0 = 0.

# Driving force: Amplitude / drive strength (gamma)  / Frequency
gamma = 0.9
fD = 1.0   
wD = 2*pi*fD
TD = 1/fD

# Time span: t1 to t2  / restircted time span t[NS]
N = 5999          
t1 = 0.0
t2 = 20
NS = -1000    # Frequency spectrum
ns = -1000        # Start time for restricted phase space plot
t = np.linspace(t1,t2,N)
dt = t[1] - t[0]

# Natural frequency and period
w0 = 1.5*wD
f0 = w0/(2*pi) 
T0 = 1/f0

# Damping constant  beta --> b
b =  w0/4
# b = w0/8
b = 3
# Pendulum length
g = 9.8
L = g/w0**2

 
#%%  SOLVE ODE
# Solve ODE for x,y    x = theta   y = omega 
def lorenz(t, state):    
    x, y = state
    dx = y
    dy = -w0**2*sin(x) - 2*b*y + gamma*w0**2*cos(wD*t)
    return [dx, dy]  

u0 = [theta0, omega0]
sol = odeint(lorenz, u0, t, tfirst=True)
theta = sol[:,0]     # angular displacement [rad/pi]
omega = sol[:,1]        # angular velocity  [rad/s]      
 

#%%
# Wrapping theta
thetaW = np.arctan2( np.sin(theta), np.cos(theta))
# theta  -pi to +pi
thetaP = np.zeros(N)
for c in range(N):
     thetaP[c] = theta[c]
     while thetaP[c] > pi:
           thetaP[c] = thetaP[c] - 2*pi
     while thetaP[c] < -pi:
          thetaP[c] = thetaP[c] + 2*pi
        
x =  L*sin(theta)    # Horizontal displacement  x


#%%   Find peaks in x vs t plot
flagP = 1.0; Tpeaks = 0.0; fpeaks = 0.0; numT = 0.0
#peaks, _ = find_peaks(theta)
peaks, _ = find_peaks(x)
z = len(peaks)
if z < 1:
   flagP = 0.0
if flagP > 0.1:   
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
Fmax = 6
Fmin = 0
nF = 2999
F = np.linspace(Fmin,Fmax,nF)
HR = np.zeros(nF); HI = HR; H = HR; H = HR+HR*1j
for c in range(nF):
     g = omega[NS:-1]*np.exp(p*F[c]*t[NS:-1])
     #g = x[NS:-1]*np.exp(p*F[c]*t[NS:-1])
     H[c] = simpson1d(g,t1,t2)

psd = 2*np.conj(H)*H   
psd = np.real(psd/max(psd))   

Find = np.where(np.real(psd)>= 0.98)
II = Find[-1]
II = II[-1]
Fpeak = abs(F[II])
Tpeak = 1/Fpeak

  
#%% CONSOLE OUTPUT
print(' ')
p = theta[0]/pi
print('Model Parameters')
print('  theta(0)/pi = %2.3f ' % p, '  omega(0) = %2.4f rad/s' % omega0 \
     ,'  Damping: b = %2.3f' % b )
print('Driving force')
print('  gamma = %2.3f' % gamma, '  TD = %2.3f s' % TD \
      ,'  fD = %2.3f Hz' % fD, '   wD = %2.3f rad/s' %wD)    
print('Time span')
print('  time steps = %2.0f' % N, '  tMax = %2.3f s' % t2) 
print('Free vibration')
print('               T0 = %2.3f s' % T0, '  f0 = %2.3f Hz' % f0,  '   w0 = %2.3f rad/s' % w0 )

print('Results')
print('Peaks t vs x:  Tpeak = %2.3f s'  % Tpeaks, '  fPeaks = %2.3f Hz' %fpeaks)
print('psd:           TPeak = %2.3f s' % Tpeak,   '  fPeak = %2.3f Hz' % Fpeak)

      
#%%
# GRAPHICS
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12


#%% Fig 0   Time evolution theta
plt.rcParams["figure.figsize"] = (4,3)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.91, bottom = 0.18, left = 0.180,\
                    right = 0.96, hspace = 0.36,wspace=0.50)

R = 0;   # t Vs theta 
axes.set_ylabel(r'$\theta$ / $\pi$',color= 'black',fontsize = 12)
axes.set_xlabel('$t$  [ s ]',color = 'black',fontsize = 12)
axes.set_title('$\gamma = $ %2.4f ' % gamma )
     
#axes.set_xlim([20, 40])
#axes[R].set_ylim([0,1])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_yticks(np.arange(0,1.1,0.2))
axes.xaxis.grid()
axes.yaxis.grid()
xP = t; yP = theta/pi
axes.plot(xP, yP, 'b')

plt.savefig('a0.png') 

#%% Fig 1  Time evolution      
plt.rcParams["figure.figsize"] = (8,7)
fig1, axes = plt.subplots(nrows=3, ncols=1)
fig1.subplots_adjust(top = 0.98, bottom = 0.12, left = 0.180,\
                    right = 0.95, hspace = 0.36,wspace=0.50)

R = 0;   # t Vs theta 
axes[R].set_ylabel(r'$\theta$ / $\pi$',color= 'black',fontsize = 12)
axes[R].set_xlabel('$t$  [ s ]',color = 'black',fontsize = 12)
#axes[R].set_xlim([0, 2])
#axes[R].set_ylim([0,1])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_yticks(np.arange(0,1.1,0.2))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
xP = t; yP = theta/pi
axes[R].plot(xP, yP, 'blue')

R = 1;   # t vs omega 
axes[R].set_ylabel(r'$\omega$  [ rad.$s^{-1}$ ]',color= 'black',fontsize = 12)
axes[R].set_xlabel('$t$  [ s ]',color = 'black',fontsize = 12)
#axes[R].set_xlim([0, 2])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_ylim([-20, 80])
#axes[R].set_yticks(np.arange(-20,81,20))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
xP = t; yP = omega
axes[R].plot(xP, yP, 'blue')

R = 2;   # t vs x 
axes[R].set_ylabel(r'$x/L$  ',color= 'black',fontsize = 12)
axes[R].set_xlabel('$t$  [ s ]',color = 'black',fontsize = 12)
#axes[R].set_xlim([0, 2])
#axes[R].set_xticks(np.arange(0,2.1,0.2))
#axes[R].set_ylim([-20, 80])
#axes[R].set_yticks(np.arange(-20,81,20))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
xP = t; yP = x/L
axes[R].plot(xP, yP, 'blue')
if gamma > 0:
   xP = t; yP = cos(wD*t)
#  axes[R].plot(xP, yP, 'r',lw = 1)

plt.savefig('a1.png') 

#%% Fig 2  Phase space plots
            

plt.rcParams["figure.figsize"] = (5,5)
fig1, axes = plt.subplots(nrows=2, ncols=1)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.180,\
                    right = 0.95, hspace = 0.36,wspace=0.40)

R = 0;   # theta vs omega
axes[R].set_ylabel(r'$\omega$ [rad.$s{^1}$]',color= 'black',fontsize = 12)
# axes[R].set_xlabel(r'$\theta$ / $\pi$',color= 'black',fontsize = 12)
# axes[R].set_xlim([0, 2])
# axes[R].set_ylim([0,1])
# axes[R].set_xticks(np.arange(0,2.1,0.2))
# axes[R].set_yticks(np.arange(0,1.1,0.2))
axes[R].xaxis.grid(); axes[R].yaxis.grid()
xP = theta/pi; yP = omega
axes[R].plot(xP, yP, 'blue')
xP = theta[0]/pi; yP = omega[0]
axes[R].plot(xP, yP, 'go', ms = 8)
xP = theta[-1]/pi; yP = omega[-1]
axes[R].plot(xP, yP, 'ro', ms = 8)

R = 1;   # theta vs omega restricted
#NS = -1000   # for tStart index
axes[R].set_ylabel(r'$\omega$ [rad.$s{^1}$]',color= 'black',fontsize = 12)
axes[R].set_xlabel(r'$\theta$ / $\pi$',color= 'black',fontsize = 12)
axes[R].set_title('t$_S$ = %2.1f s' % t[ns] + '   t$_F$ = %2.1f s' % t[-1]  \
         , fontsize = 12)
# axes[R].set_xlim([0, 2])
# axes[R].set_ylim([0,1])
# axes[R].set_xticks(np.arange(0,2.1,0.2))
# axes[R].set_yticks(np.arange(0,1.1,0.2))
axes[R].xaxis.grid(); axes[R].yaxis.grid()


xP = theta[ns:-1]/pi; yP = omega[ns:-1]
#xP = thetaP[NS:-1]/pi; yP = omega[NS:-1]
axes[R].plot(xP, yP, 'blue')

xP = theta[0]/pi; yP = omega[0]
axes[R].plot(xP, yP, 'go', ms = 8)
xP = theta[-1]/pi; yP = omega[-1]
axes[R].plot(xP, yP, 'ro', ms = 8)

plt.savefig('a2.png') 

#%% Fig 3  Frequency spectrum 
plt.rcParams["figure.figsize"] = (5,5)
fig1, axes = plt.subplots(nrows=2, ncols=1)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.180,\
                    right = 0.95, hspace = 0.36,wspace=0.40)

R = 0;   # f vs psd 
axes[R].set_title('t$_S$ = %2.1f s' % t[NS] + '   t$_F$ = %2.1f s' % t[-1]  \
         , fontsize = 12)
axes[R].set_ylabel('$psd$',color= 'black',fontsize = 12)
axes[R].set_xlabel('$f$  [ Hz ]',color = 'black',fontsize = 12)
axes[R].set_xlim([0, Fmax])
axes[R].set_ylim([0,1])
axes[R].set_xticks(np.arange(0,6.1,0.5))
axes[R].set_yticks(np.arange(0,1.1,0.2))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
xP = F[F>0]; yP = psd[F>0]
axes[R].plot(xP, yP, 'blue')
xP = [f0,f0]; yP = [0,1]
axes[R].plot(xP, yP, 'r')
xP = [fD,fD]; yP = [0,1]
axes[R].plot(xP, yP, 'm')

R = 1;   # f vs log(psd) 
axes[R].set_ylabel('$log(psd)$',color= 'black',fontsize = 12)
axes[R].set_xlabel('$f$  [ Hz ]',color = 'black',fontsize = 12)
axes[R].set_xlim([0, 6])
axes[R].set_xticks(np.arange(0,6.1,0.5))
#axes[R].set_ylim([-20, 80])
#axes[R].set_yticks(np.arange(-20,81,20))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
xP = F[F>0]; yP = np.log(psd[F>0])
axes[R].plot(xP, yP, 'blue')
xP = [f0,f0]; yP = [-20,0]
axes[R].plot(xP, yP, 'r')
xP = [fD,fD]; yP = [-20,0]
axes[R].plot(xP, yP, 'm')

plt.savefig('a3.png')    

#%%  Find and plot peaks

RP = np.arange(3000,5998,1)
#y = theta[RP]/pi
y = omega[RP]

tP = t[RP]
# Findpeaks
L = len(y)
a1 = y[0]; a2 = y[1]
flag = 2
indexMax = np.zeros(20)

if a2 > a1:
    flag = 1
v = 0
for c in range(L-1):
    y1 = y[c]     
    y2 = y[c+1]

    if flag == 1 and y2 > y1:
       c = c+1
    if (flag == 1 and y2 < y1):
       indexMax[v] = c
       v = v + 1
       c = c+1
    if y2 <= y1:
      flag = 0
    if y2 > y1:
      flag = 1

ind = indexMax[indexMax > 0]
ind = ind.astype(int)


font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize = (3, 2))

fig.subplots_adjust(top = 0.94, bottom = 0.27, left = 0.22,\
                     right = 0.96, hspace = 0.2,wspace=0.2)
   
xP = tP; yP = y
plt.plot(xP,yP,linewidth=2,color='b')
xP = tP[ind]
yP = y[ind]
plt.plot(xP,yP,'ro',ms = 6)
plt.ylim([0,max(omega)])
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
plt.grid('visible')
plt.xlabel(r'$t$   [ s ]', fontdict = font1)
plt.ylabel(r'$\omega$ [rad.$s{^1}$]',color= 'black',fontsize = 12)
#plt.ylabel(r'$\theta$ / $\pi$',color= 'black',fontsize = 12)
plt.savefig('a4.png')    

for c in range(len(yP)-1):
    print(r'   %2.3f ' %xP[c], '    %2.3f  ' %yP[c])

#%%
tEnd = time.time()
tE = tEnd - tStart


print('\nExecution time =  %2.3f s' % tE )


