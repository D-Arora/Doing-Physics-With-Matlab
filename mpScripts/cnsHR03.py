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

# FUNCTIONS  ==============================================================
# Solve ODE for x,y,z
def lorenz(t, state, k):    
    x, y, z = state
    
    dx = y - x**3 + 3*x**2 - z + Iext 
    dy = 1 - 5*x**2 - y
    dz = r * ( s*(x - x0) - z )
    
    return [dx, dy, dz]

#%%   
# CONSTANTS  ===============================================================

Iext = 1.5

r = 0.005

s = 4



k = np.zeros(2)


# SETUP  ============================================================


# >>>>> Time interval t from t1 to t2 / plot limits
t1 = 0
t2 = 2500
N = 9999
R1 = N-100
R2 = N-1
# Initial values x0, y0, z0
x0 = -1.6
u0x = x0
u0y = 0
u0z = 0
u0 = [u0x, u0y, u0z]

# Time span and limits
         
tSpan = np.linspace(t1,t2,N)
t = tSpan


# SOLVE ODEs  ============================================================
sol = odeint(lorenz, u0, tSpan, args = (k,), tfirst=True)
x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

e = 1e-3

u0x = u0x + e
u0y = u0y + e
u0z = u0z + e
u0 = [u0x, u0y, u0z]

sol = odeint(lorenz, u0, tSpan, args = (k,), tfirst=True)
X = sol[:,0]
Y = sol[:,1]
Z = sol[:,2]

dx = abs(X[R1:R2]-x[R1:R2])



S = sum(dx)/e
L = np.log(S)/len(dx)
print(L)


dfx = -3*x[R1:R2]**2 + 6*x[R1:R2]
tR = t[R1:R2]
yR = np.log(abs(dfx))

plt.plot(tR,yR)
plt.ylim([-5,5])

plt.plot(tR,np.log(abs(dx/e)))

# Fig 1  x vs t  ---------------------------------------------------------
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize = (5, 4))
fig.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.15,\
                    right = 0.92, hspace = 0.60,wspace=0.5)

plt.plot(t[R1:R2],x[R1:R2]+1,linewidth=2,color='b')
plt.plot(t[R1:R2],X[R1:R2],linewidth=2,color='r')
plt.plot(t[R1:R2],dx,linewidth=2,color='m')
#plt.xlim([R1,R2])
plt.ylim([-2,2])
#yR = np.arange(0,250,100) 
#plt.yticks(yR)
plt.grid('visible')
plt.ylabel('x', fontdict = font1)
plt.xlabel('t ', fontdict = font1)

# Fig 1  dx vs t  ---------------------------------------------------------
# font1 = {'family':'Tahoma','color':'black','size':12}
# plt.rcParams['font.family'] = ['Tahoma']
# plt.rcParams['font.size'] = 12

# fig = plt.figure(figsize = (5, 4))
# fig.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.15,\
#                     right = 0.92, hspace = 0.60,wspace=0.5)


# plt.plot(t[R1:R2]-t[R1],np.log(dx/e),linewidth=2,color='m')
# #plt.xlim([R1,R2])
# #plt.ylim([-2,2])
# #yR = np.arange(0,250,100) 
# #plt.yticks(yR)
# plt.grid('visible')
# plt.ylabel('x', fontdict = font1)
# plt.xlabel('t ', fontdict = font1)

