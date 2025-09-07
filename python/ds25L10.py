# -*- coding: utf-8 -*-
"""

Bifurcations in a Model for Insect Outbreak - Dynamical Systems 

# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/ds25L10.pdf

"""

# LIBRARIES  ================================================================
import numpy as np
from numpy import pi, sin, cos, linspace
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')


#%%  SOLVE ODE for x   where x is the popiulation N
def lorenz(t, state):    
    x = state
    dx = r*x*(1 - x/K) - B*x**2/(A**2 + x**2)
    return dx  

#%%
r,K = 1.5, 1.0
A,B = 0.118, 0.3               # A = 0.01   B = 0.08   
num, tMax = 999,30
t = linspace(0,tMax,num)



#%%
# Fig 1: Membrane potential x vs time t
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('N')
ax.set_title('A = %0.3f' %A + '  B = %0.3f' %B)
ax.set_ylim([0,1.1])
ax.grid()

x0 = np.array([0.1,0.2,0.3,1.0])

for c in range(len(x0)):
    u0 = x0[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    xS = sol[:,0] 
    ax.plot(t,xS,'b',lw = 2)

fig1.tight_layout()

#%% Fig 2: x vs xDot
x = linspace(0,0.8*K,num)
xDot = r*x*(1 - x/K) - B*x**2/(A**2 + x**2)

# Find zeros in Ndot
Q = np.zeros(5); p = 0
for c in range(num-2):
    q = xDot[c]*xDot[c+1]
    if q <= 0:
       Q[p] = c
       p = int(p+1)     
QI = Q.astype(int)
xZ = x[QI]          # Zeros for xDot  
print(xZ)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('N',fontsize = 14)
ax.set_ylabel('N$_{dot}$',fontsize = 14)
ax.set_title('A = %0.3f' %A + '  B = %0.3f' %B)
ax.grid()
ax.plot(x,xDot,'k',lw = 2)
ax.plot([0,max(x)],[0,0],'m',lw = 1.6)
#ax.plot(xZ[1],0,'bo',ms=7)
#ax.arrow(0.2, 0, 0.2, 0,  lw = 3,head_width=0.03, head_length=0.05, fc='k', ec='k')
#ax.arrow(1.4, 0, -0.2, 0,  lw = 3,head_width=0.03, head_length=0.05, fc='k', ec='k')
fig2.tight_layout()


#%%  Fig 3: slope plot  x vs t
N = 15
tQ = linspace(0,30,N)
x = linspace(0,1,N)

f = r*x*(1 - x/K) - B*x**2/(A**2 + x**2)

T,X = np.meshgrid(tQ,x)

dX = np.ones([N,N])
F =  r*X*(1 - X/K) - B*X**2/(A**2 + X**2)

dY = F/(np.sqrt(dX**2 + F**2))
dX = dX/(np.sqrt(dX**2 + F**2))

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('t')
ax.set_ylabel('N')
ax.set_title('A = %0.3f' %A + '  B = %0.3f' %B)
ax.set_ylim([0,1])
ax.set_xlim([0,max(tQ)])

ax.quiver(T,X,dX,6*dY,color = 'b')

for c in range(len(x0)):
    u0 = x0[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    xS = sol[:,0] 
    ax.plot(t,xS,'r',lw = 2)

fig3.tight_layout()  

#%%  Fig 4: SLOPE FIELD
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('t')
ax.set_ylabel('N')
ax.set_ylim([0,1])
ax.set_xlim([0,max(tQ)])

ax.streamplot(T,X,dX,dY,color = 'b')
for c in range(len(x0)):
    u0 = x0[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    xS = sol[:,0] 
    ax.plot(t,xS,'r',lw = 2)

fig4.tight_layout()  

#%% Fig. 5:  PREDATION CURVE  P
p = linspace(0,2,999)
P = B*p**2/(A**2+p**2)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('N',fontsize = 12)
ax.set_ylabel('N$_{dot}$',fontsize = 12)
ax.set_title('A = %0.3f' %A + '  B = %0.3f' %B)
ax.grid()
ax.plot(p,P,'k',lw = 2)
fig5.tight_layout()


#%% Fig 6:  x vs y1 and y2
# A = 0.1; B = 0.3
x = linspace(0,1.1,599)
y1 = r*(1-x/K)
y2 = B*x/(A**2 + x**2)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig6, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('N',fontsize = 14)
ax.set_ylabel('y ',fontsize = 14)
ax.set_title('A = %0.3f' %A + '  B = %0.3f' %B)
ax.grid()
ax.plot(x,y1,'r',lw = 2,label = 'y$_1$')
ax.plot(x,y2,'b',lw = 2,label = 'y$_2$')
ax.legend()

fig6.tight_layout()


#%%         
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')
fig6.savefig('a6.png')

'''

