# cs_005.py
# Ian Cooper
# Jan 2024
# COMPLEX SYSTEMS
#	LINEAR DISCRETE-TIME DYNAMICAL SYSTEMS:
#       ASYMPTOTIC BEHAVIOUR – EIGENVALUES AND EIGENVECTORS
  

# # Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_005.pdf

# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.signal import chirp, find_peaks, peak_widths
import time
from numpy.linalg import eig


# %% Fibonacci sequence
# >>> Input number of time steps 
nT = 99
xR = np.zeros(nT)
# >>> Input initial conditions
xR[0] = 0; xR[1] = 1
# Recursion relationship
for k in range(nT-2):
    xR[k+2] = xR[k+1] + xR[k]

tR = np.arange(0,nT,1)

#%%  x(t+1) = A x(t)
# >>> Transformation matrix
A = np.array([[1,1],[1,0]])
# >>> Initial conditions
X0 = np.array([0,1])
# Eigenvalues and eigenvectors
e, f = eig(A)
# Eigenvectors 
v0 = f[:,0]
v1 = f[:,1]
# Solve simutaneous equations
M = np.array([ [v0[0],v1[0]],[v0[1],v1[1]] ])
c = np.linalg.solve(M,X0)

# Time evolution
N = 99
x = np.zeros(N); y = np.zeros(N)
t = np.arange(N)
x[0] = 1; x[1] = 2
L0 = max(e)
L = e/L0

for k in range(N-1):
  # x[k+1] = x[k] + c[0]*L[0]**k*v0[0] + c[1]*L[1]**k*v1[0]
  #  y[k+1] = y[k] + c[0]*L[0]**k*v0[1] + c[1]*L[1]**k*v1[1]
  x[k] = c[0]*e[0]**k*v0[0] + c[1]*e[1]**k*v1[0]
  y[k] = c[0]*e[0]**k*v0[1] + c[1]*e[1]**k*v1[1]
#x = L0*x
#y = L0*x

# Time evolution using largest eigenvalue
nE = 19
xE = np.zeros(nE)
tE = np.arange(0,nE)
for k in range(nE-1):
    xE[k+1] = L0*xE[k]

    
    
#%% CONSOLE OUTPUT
print(' ')
print('Fabonacci sequence ')
print('Initial conditions: x0 = %2.0f ' % abs(x[0]), 'x1 = %2.0f' % x[1])
print(xR[0:9])
print('Matrix A')
print(A)
print('Eigenvalues e')
print('e0 = %5.5f ' % e[0], 'e1 = %5.5f' % e[1])
print('Eigenfunctions v:  v0   v1')
print(v0,v1)
print('c coefficents')
print('c0 = %5.5f ' % c[0], 'c1 = %5.5f' % c[1])
      
#%%
# GRAPHICS
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

def graph():
    fig = plt.figure(figsize = (4, 3))
    fig.subplots_adjust(top = 0.94, bottom = 0.20, left = 0.15,\
                        right = 0.92, hspace = 0.60,wspace=0.5)
    return

# Fig 1  x vs t  ---------------------------------------------------------
graph()

plt.plot(tR,xR,'bo')
plt.plot(t,x,linewidth=2,color='r')
plt.plot(tE,xE,linewidth=2,color='m')
#yR = np.arange(0,250,100) 
#plt.yticks(yR)
plt.grid('visible')
plt.ylabel('x', fontdict = font1)
plt.xlabel('t ', fontdict = font1)
#plt.savefig('cns001.png')


# fig. 2.  Eigenfunctions v0 and v1
fig = plt.figure(figsize = (4, 4))
ax = fig.add_subplot()

 
# square plot
ax.set_aspect('equal', adjustable='box')
fig.subplots_adjust(top = 0.96, bottom = 0.17, left = 0.18,\
                    right = 0.96, hspace = 0.60,wspace=0.5)
xP = [0,v0[0]]; yP = [0,v0[1]]
plt.plot(xP,yP,linewidth=2,color='b',label = v0)
xP = [0,v1[0]]; yP = [0,v1[1]]
plt.plot(xP,yP,linewidth=2,color='r', label = v1)

plt.xlim([-1,1])
xR = np.arange(-1,1.1,0.5) 
plt.xticks(xR)
plt.ylim([-1,1])
yR = np.arange(-1,1.1,0.5) 
plt.grid('visible')
plt.legend()

plt.savefig('cs001.png')



# A = np.array([ [1,2], [3,4] ])
# R = A@A@A@A@A@A

# print(A)
# print(R)             
             