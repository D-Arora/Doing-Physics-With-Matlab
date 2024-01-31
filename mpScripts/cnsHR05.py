# -*- coding: utf-8 -*-
# cnsHR01.py
# 22 Dec 2023
# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# HINDMARSH - ROSE
# SPIKING and BURSTING NEURON


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9871131/

# https://scicomp.stackexchange.com/questions/36013/numerical-computation-of-lyapunov-exponent


# LIBRARIES  ==============================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import eig
from numpy import real, imag
from scipy.signal import chirp, find_peaks, peak_widths
import time
# FUNCTIONS  ==============================================================
# Solve ODE for x,y,z
def lorenz(t, state, k):    
    x, y, z = state
    
    dx = y - x**3 + 3*x**2 - z + Iext 
    dy = 1 - 5*x**2 - y
    dz = r * ( s*(x - x0) - z )
    
    return [dx, dy, dz]

#%% 
time1 = time.time()
# Fig 1  ISI vs I  ---------------------------------------------------------
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize = (5, 4))
fig.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.15,\
                    right = 0.92, hspace = 0.60,wspace=0.5)

#yP = ISI
#xP = np.zeros(L-1)+Iext
#plt.plot(xP,yP,'o',ms = 2)
plt.xlim([1,4])
plt.ylim([0,200])
#yR = np.arange(0,250,100) 
#plt.yticks(yR)
plt.grid('visible')
plt.ylabel('ISI', fontdict = font1)
plt.xlabel('I ', fontdict = font1)
  
# CONSTANTS  ===============================================================



r = 0.005

s = 4

k = np.zeros(2)


# SETUP  ============================================================
# >>>>> External current


# >>>>> Time interval t from t1 to t2 / plot limits
t1 = 0
t2 = 1000
R1 = 000; R2 = t2

# Initial values x0, y0, z0
x0 = -1.6
u0x = 0.1
u0y = 1
u0z = 0.2
u0 = [u0x, u0y, u0z]

# Time span and limits
N = 9999            
tSpan = np.linspace(t1,t2,N)
t = tSpan
T1 = 6666
T2 = N

Imin = 1; Imax = 4; num = 999     # 9999
I = np.linspace(Imin,Imax,num)

for c in range(num):
# SOLVE ODEs  ============================================================
   Iext = I[c]
   sol = odeint(lorenz, u0, tSpan, args = (k,), tfirst=True)
   x = sol[:,0]
   y = sol[:,1]
   z = sol[:,2]

# # Find peaks
   peaks, _ = find_peaks(x)
   tPeaks = t[peaks]
   tP = tPeaks[tPeaks>300]
   L = len(tP)
   ISI = np.zeros(L-1)

   for c in range(L-1):
      ISI[c] = tP[c+1]-tP[c]
    
   yP = ISI
   xP = np.zeros(L-1)+Iext
   plt.plot(xP,yP,'ob',ms = 0.25)
            
            
   
time2 = time.time()

print((time2-time1)/60)

# CONSOLE OUTPUT  ===========================================================

# print('Model parameters')
# print('r = %5.4f ' % r, 's = %5.4f' % s, 'Iext = %5.4f ' % Iext)
# print('Initial values')
# print('x0 = %5.4f ' % u0x, 'y0 = %5.4f' % u0y, 'z0 = %5.4f' % u0z)
# print('  ')
#print('ISI')
#print(ISI)

# for c in range(3):
#     flag1 = 0
#     print('xC = %2.3f ' % real(xC[c]), '%2.3f' % imag(xC[c]),'j', \
#        '   yC = %2.3f ' % real(yC[c]), '%2.3f' % imag(yC[c]),'j', \
#        '   zC = %2.3f ' % real(zC[c]), '%2.3f' % imag(zC[c]),'j')
#     print('eV = %2.3f ' % real(eV[c]), '%2.3f' % imag(eV[c]),'j') 
#     if imag(xC[c]) == 0 and imag(yC[c]) == 0 and imag(zC[c]) == 0 and real(eV[c]) < 0:
#        print('Stable fixed point')    
#     else:
#        print('NOT stable fixed point')   
#     print('  ')


#%%

#%%           
# GRAPHICS ================================================================
#plt.close('all')



# Fig 2  x vs t  ---------------------------------------------------------
# font1 = {'family':'Tahoma','color':'black','size':12}
# plt.rcParams['font.family'] = ['Tahoma']
# plt.rcParams['font.size'] = 12

# fig = plt.figure(figsize = (5, 3))
# fig.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.15,\
#                     right = 0.92, hspace = 0.60,wspace=0.5)

# plt.plot(t,x,linewidth=2,color='b')
# plt.xlim([R1,R2])
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
# plt.grid('visible')
# plt.ylabel('x', fontdict = font1)
# plt.xlabel('t ', fontdict = font1)
#plt.savefig('cns001.png')

# # Fig 2  y vs t  ---------------------------------------------------------

# fig = plt.figure(figsize = (5, 3))
# fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.2,\
#                     right = 0.92, hspace = 0.60,wspace=0.5)

# plt.plot(t,y,linewidth=2,color='b')
# plt.xlim([R1,R2])
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
# plt.grid('visible')
# plt.ylabel('y', fontdict = font1)
# plt.xlabel('t ', fontdict = font1)
# plt.savefig('cns002.png')

# # Fig 3  z vs t  ----------------------------------------------------------
# fig = plt.figure(figsize = (5, 3))
# fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.15,\
#                     right = 0.92, hspace = 0.60,wspace=0.5)

# plt.plot(t,z,linewidth=2,color='b')
# #plt.xlim([2500,3500])
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
# plt.grid('visible')
# plt.ylabel('z', fontdict = font1)
# plt.xlabel('t ', fontdict = font1)
# plt.savefig('cns003.png')

# #%% Fig. 4    PHASE PORTRAIT PLOT  x vs y
# fig = plt.figure(figsize = (5, 3))
# fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.20,\
#                     right = 0.96, hspace = 0.60,wspace=0.5)

# # Stream plot or quiver plot
# X = np.arange(-2.5, 3.5, 1)
# Y = np.arange(-25, 20, 4)
 
# XX, YY = np.meshgrid(X, Y)
# ZZ = z[N-1]
# II = Iext
# u = YY - XX**3 + 3*XX**2 - ZZ + II
# v = 1 - 5*XX**2 - YY

# us = u/(u**2 + v**2)**0.5
# vs = v/(u**2 + v**2)**0.5
# #plt.quiver(XX, YY, us, vs)
# plt.streamplot(XX, YY, us, vs)    
 
# # Plot nullclines   
# plt.plot(xN,yX,linewidth=2,color='r')
# plt.plot(xN,yY,linewidth=2,color='m')

# # Plot phase portrait
# plt.plot(x[T1:T2],y[T1:T2],linewidth=3,color='b')
# #plt.ylim([0,250])
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
# plt.grid('visible')
# plt.ylabel('y', fontdict = font1)
# plt.xlabel('x ', fontdict = font1)

# # Plot stable fixed point (xc, yC)
# if imag(xC[0]) == 0:
#     plt.plot(xC[0],yC[0],'ro',ms = 8)
# if imag(xC[1]) == 0:
#     plt.plot(xC[1],yC[1],'ro',ms=8)
# if imag(xC[2]) == 0:
#     plt.plot(xC[2],yC[2],'ro',ms=8)    

# plt.savefig('cns004.png')    

# # Fig 5  x vs z  ----------------------------------------------------------
# fig = plt.figure(figsize = (5, 3))
# fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.15,\
#                     right = 0.92, hspace = 0.60,wspace=0.5)

# plt.plot(x,z,linewidth=2,color='b')
# plt.plot(x[N-1],z[N-1],'ro',ms = 8)
# #plt.xlim([2500,3500])
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
# plt.grid('visible')
# plt.ylabel('z', fontdict = font1)
# plt.xlabel('x ', fontdict = font1)
# plt.savefig('cns005.png')

# #%% Fig. 6  PLOT3D   x,y,z vs t  ------------------------------------------
# fig = plt.figure(figsize = (5, 3))
# fig.subplots_adjust(top=0.98, bottom = 0.12, left = 0.01,\
#                     right = 0.96, hspace = 0.60,wspace=0.5)

# ax = plt.axes(projection='3d')

# ax.plot3D(x[T1:T2],y[T1:T2],z[T1:T2])
# ax.plot3D(x[N-1],y[N-1],z[N-1],'ro')
# plt.xlabel('x')
# plt.ylabel(' y')
# ax.set_zlabel('  z');
# plt.savefig('cns006.png')
