# -*- coding: utf-8 -*-
"""
cs101.py
Aug 25

DYNAMICAL SYSTEMS:
    Transcritical Bifurcations: subcritical

# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L7.pdf

"""


# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace, exp

plt.close('all')

tStart = time.time()

#%%   FUNCTIONS
def funct(a,b,x):
    xDot = (1-x)*x - a*(1-exp(-b*x))
    return xDot

#%%  SOLVE ODE
# Solve ODE for x
def lorenz(t, state):    
    x = state
    dx = (1-x)*x - a*(1-exp(-b*x))
    return dx  

#%%  Find zeros of xDot
def ZEROS():
    p = 0; Q = np.zeros(2)
    for c in range(num-2):
        q = xDot[c]*xDot[c+1]
        if q <= 0:
           Q[p] = c
           p = int(p+1)     
    Q = Q.astype(int)
    Z = x[Q]
    return Z



#%%  Fig 1: Plot x vs xDot
num = 599
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,3)
fig1, ax = plt.subplots(1,3)

C = 0; a = 0.9; b = 1
x1, x2 = -0.1, 0.3
x = linspace(x1,x2,num)
xDot = funct(a,b,x)
xZ = ZEROS()                # zeros xDot

ax[C].set_xlabel('x',color = 'black')
ax[C].set_ylabel('f(x) = $x_{dot}$',color= 'black')   
ax[C].grid()
ax[C].set_title('a = %0.1f' %a + '   b = %0.1f' %b, fontsize = 11)
ax[C].plot(x,xDot,'k',lw = 2)
ax[C].plot(xZ[0],0,'ro',ms = 6)
ax[C].plot(xZ[1],0,'bo',ms = 6)

C = 1; a = 1.0; b = 1
x1, x2 = -0.2,0.2
x = linspace(x1,x2,599)
xDot = funct(a,b,x)
ax[C].set_xlabel('x',color = 'black')
ax[C].xaxis.grid()
ax[C].yaxis.grid()
ax[C].set_title('a = %0.1f' %a + '   b = %0.1f' %b, fontsize = 11)
ax[C].plot(x,xDot,'k',lw = 2)
ax[C].plot(0.01,0,'bo',ms = 6)
ax[C].plot(-0.01,0,'ro',ms = 6)

C = 2; a = 1.1; b = 1
x1, x2 = -0.3,0.1
x = linspace(x1,x2,599)
xDot = funct(a,b,x)
xZ = ZEROS()
ax[C].set_xlabel('x',color = 'black')
ax[C].xaxis.grid()
ax[C].yaxis.grid()
ax[C].set_title('a = %0.1f' %a + '   b = %0.1f' %b, fontsize = 11)
xP = x; yP = xDot; 
ax[C].plot(xP,yP,'k',lw = 2)
ax[C].plot(xZ[0],0,'ro',ms = 6)
ax[C].plot(xZ[1],0,'bo',ms = 6)
fig1.tight_layout()

# xxx
# #%%  Fig 2:
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (7,3)
# fig2, ax = plt.subplots(1,2)

# r = 0; xC = 0
# C = 0
# x0 = 0.5
# tMax = 60 
# t = linspace(0,tMax,999)
# sol = odeint(lorenz, x0, t, tfirst=True)
# xS = sol[:,0] 
# ax[C].set_xlabel('t',color = 'black')
# ax[C].set_ylabel('x',color= 'black')
# ax[C].set_title(r'x(0) = %2.2f' %x0, fontsize = 12)
# ax[C].grid()
# ax[C].plot(t,xS,'b',lw = 2)    

# C = 1
# x0 = -0.1
# tMax = 9.6 
# t = linspace(0,tMax,999)
# sol = odeint(lorenz, x0, t, tfirst=True)
# xS = sol[:,0] 
# ax[C].set_xlabel('t',color = 'black')
# ax[C].set_ylabel('x',color= 'black')
# ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
# ax[C].grid()
# ax[C].plot(t,xS,'r',lw = 2) 
# fig2.suptitle('r = %0.0f' %r  + '    x$_e$ = 0')
# fig2.tight_layout() 

# #%%  Fig 3:
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (7,3)
# fig3, ax = plt.subplots(1,3)

# r = -10
# xC = np.zeros(2); xC[0] = 0; xC[1] = r
# C = 0
# x0 = -10.5
# tMax = 0.3
# t = linspace(0,tMax,999)
# sol = odeint(lorenz, x0, t, tfirst=True)
# xS = sol[:,0] 
# ax[C].set_xlabel('t',color = 'black')
# ax[C].set_ylabel('x',color= 'black')
# ax[C].set_title(r'x(0) = %2.2f' %x0, fontsize = 12)
# ax[C].grid()
# ax[C].plot(t,xS,'r',lw = 2)    

# C = 1
# x0 = -9.5
# tMax = 2.1
# t = linspace(0,tMax,999)
# sol = odeint(lorenz, x0, t, tfirst=True)
# xS = sol[:,0] 
# ax[C].set_xlabel('t',color = 'black')
# ax[C].set_ylabel('x',color= 'black')
# ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
# ax[C].grid()
# ax[C].plot(t,xS,'b',lw = 2) 

# C = 2
# x0 = 5
# tMax = 2.1
# t = linspace(0,tMax,999)
# sol = odeint(lorenz, x0, t, tfirst=True)
# xS = sol[:,0] 
# ax[C].set_xlabel('t',color = 'black')
# ax[C].set_ylabel('x',color= 'black')
# ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
# ax[C].grid()
# ax[C].plot(t,xS,'b',lw = 2) 

# fig3.suptitle('r = %0.0f' %r +'   x$_e$ = %0.0f' %xC[1]     + '   x$_e$ = %0.0f' %xC[0])
# fig3.tight_layout() 


# #%%  Fig 4: r =+10  t vs x
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (7,3)
# fig4, ax = plt.subplots(1,3)

# r = 10
# xC = np.zeros(2); xC[0] = 0; xC[1] = r
# C = 0
# x0 = 12
# tMax = 1
# t = linspace(0,tMax,999)
# sol = odeint(lorenz, x0, t, tfirst=True)
# xS = sol[:,0] 
# ax[C].set_xlabel('t',color = 'black')
# ax[C].set_ylabel('x',color= 'black')
# ax[C].set_title(r'x(0) = %2.2f' %x0, fontsize = 12)
# ax[C].grid()
# ax[C].plot(t,xS,'b',lw = 2)    

# C = 1
# x0 = 5
# tMax = 2.1
# t = linspace(0,tMax,999)
# sol = odeint(lorenz, x0, t, tfirst=True)
# xS = sol[:,0] 
# ax[C].set_xlabel('t',color = 'black')
# ax[C].set_ylabel('x',color= 'black')
# ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
# ax[C].grid()
# ax[C].plot(t,xS,'b',lw = 2) 

# C = 2
# x0 = -0.1
# tMax = 0.45
# t = linspace(0,tMax,999)
# sol = odeint(lorenz, x0, t, tfirst=True)
# xS = sol[:,0] 
# ax[C].set_xlabel('t',color = 'black')
# ax[C].set_ylabel('x',color= 'black')
# ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
# ax[C].grid()
# ax[C].plot(t,xS,'r',lw = 2) 

# fig4.suptitle('r = %0.0f' %r +'   x$_e$ = %0.0f' %xC[1]     + '   x$_e$ = %0.0f' %xC[0])
# fig4.tight_layout() 


# #%% r vs xE
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (4,4)
# fig5, ax = plt.subplots(1)

# ax.set_xlabel('r',color = 'black',fontsize = 12)
# ax.set_ylabel('$x_C$',color= 'black',fontsize = 14)

# ax.grid()

# ax.set_xlim([-10,10]); ax.set_ylim([-10,10])

# ax.plot([-10,0],[0,0],'b',lw = 2,label = 'Stable') 
# ax.plot([0,10],[0,0],'r',lw = 2,label = 'Ustable') 
# ax.plot([-10,0],[-10, 0],'r',lw = 2) 
# ax.plot([0,10],[0, 10],'b',lw = 2) 

# ax.legend() 
# ax.axis('equal')
# fig5.tight_layout() 

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

fig1.savefig('a1.png')

#%%
'''

fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')



'''




