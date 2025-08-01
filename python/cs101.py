# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 21:22:37 2024

@author: Owner
"""

# https://math.libretexts.org/Bookshelves/Scientific_Computing_Simulations_and_Modeling/Scientific_Computing_(Chasnov)/II%3A_Dynamical_Systems_and_Chaos/12%3A_Concepts_and_Tools

# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import time
from numpy import pi, sin, cos, linspace

plt.close('all')

tStart = time.time()

#%%   FUNCTIONS
def funct(r,x):
    xDot = r*x - x**2
    return xDot

#%%  SOLVE ODE
# Solve ODE for x
def lorenz(t, state):    
    xS = state
    dx = r*xS - xS**2
    return dx  

#%%  Fig 1: Plot x vs xDot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,3)
fig1, ax = plt.subplots(1,3)

C = 0; r = -10
x1, x2 = -20,10
x = linspace(x1,x2,599)
xDot = funct(r,x)
ax[C].set_xlabel('x',color = 'black')
ax[C].set_ylabel('f(x) = $x_{dot}$',color= 'black')   
ax[C].grid()
ax[C].set_title('r = %0.1f' %r)
ax[C].set_xlim([x1,x2])
ax[C].set_xticks([-20,-10,0,10])
ax[C].plot(x,xDot,'k',lw = 2)
ax[C].plot(0,0,'bo',ms = 6)
ax[C].plot(r,0,'ro',ms = 6)
ax[C].arrow(-12, -10, -5, 0, head_width=8, head_length=0.5, lw = 2,fc='red', ec='red')
ax[C].arrow(-8, -10, 5, 0, head_width=8, head_length=0.5, lw = 2,fc='b', ec='b')
ax[C].arrow(9, -10, -5, 0, head_width=8, head_length=0.5, lw = 2,fc='b', ec='b')

C = 1; r = 0
x1, x2 = -10,10
x = linspace(x1,x2,599)
xDot = funct(r,x)
ax[C].set_xlabel('x',color = 'black')
ax[C].set_ylabel('f(x) = $x_{dot}$',color= 'black')   
ax[C].xaxis.grid()
ax[C].yaxis.grid()
ax[C].set_title('r = %0.1f' %r)
ax[C].set_xlim([x1,x2])
ax[C].plot(x,xDot,'k',lw = 2)
ax[C].plot(0.5,0,'bo',ms = 8)
ax[C].plot(r,0,'ro',ms = 6)
ax[C].arrow(-1, -10, -5, 0, head_width=5, head_length=0.5, lw = 2,fc='red', ec='red')
ax[C].arrow(6, -10, -5, 0, head_width=5, head_length=0.5, lw = 2,fc='b', ec='b')

C = 2; r = 10
x1, x2 = -10,20
x = linspace(x1,x2,599)
xDot = funct(r,x)
ax[C].set_xlabel('x',color = 'black')
ax[C].set_ylabel('f(x) = $x_{dot}$',color= 'black')   
ax[C].xaxis.grid()
ax[C].yaxis.grid()
ax[C].set_title('r = %0.1f' %r)
ax[C].set_xlim([x1,x2])
ax[C].set_xticks([-10,0,10,20])
xP = x; yP = xDot; 
ax[C].plot(xP,yP,'k',lw = 2)
ax[C].plot(0,0,'ro',ms = 6)
ax[C].plot(r,0,'bo',ms = 6)
ax[C].arrow(-1, -10, -5, 0, head_width=8, head_length=0.5, lw = 2,fc='r', ec='r')
ax[C].arrow(3, -10, 5, 0, head_width=8, head_length=0.5, lw = 2,fc='b', ec='b')
ax[C].arrow(18, -10, -5, 0, head_width=8, head_length=0.5, lw = 2,fc='b', ec='b')

fig1.tight_layout()

#%%  Fig 2:
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,3)
fig2, ax = plt.subplots(1,2)

r = 0; xC = 0
C = 0
x0 = 0.5
tMax = 60 
t = linspace(0,tMax,999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
ax[C].set_xlabel('t',color = 'black')
ax[C].set_ylabel('x',color= 'black')
ax[C].set_title(r'x(0) = %2.2f' %x0, fontsize = 12)
ax[C].grid()
ax[C].plot(t,xS,'b',lw = 2)    

C = 1
x0 = -0.1
tMax = 9.6 
t = linspace(0,tMax,999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
ax[C].set_xlabel('t',color = 'black')
ax[C].set_ylabel('x',color= 'black')
ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
ax[C].grid()
ax[C].plot(t,xS,'r',lw = 2) 
fig2.suptitle('r = %0.0f' %r  + '    x$_e$ = 0')
fig2.tight_layout() 

#%%  Fig 3:
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,3)
fig3, ax = plt.subplots(1,3)

r = -10
xC = np.zeros(2); xC[0] = 0; xC[1] = r
C = 0
x0 = -10.5
tMax = 0.3
t = linspace(0,tMax,999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
ax[C].set_xlabel('t',color = 'black')
ax[C].set_ylabel('x',color= 'black')
ax[C].set_title(r'x(0) = %2.2f' %x0, fontsize = 12)
ax[C].grid()
ax[C].plot(t,xS,'r',lw = 2)    

C = 1
x0 = -9.5
tMax = 2.1
t = linspace(0,tMax,999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
ax[C].set_xlabel('t',color = 'black')
ax[C].set_ylabel('x',color= 'black')
ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
ax[C].grid()
ax[C].plot(t,xS,'b',lw = 2) 

C = 2
x0 = 5
tMax = 2.1
t = linspace(0,tMax,999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
ax[C].set_xlabel('t',color = 'black')
ax[C].set_ylabel('x',color= 'black')
ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
ax[C].grid()
ax[C].plot(t,xS,'b',lw = 2) 

fig3.suptitle('r = %0.0f' %r +'   x$_e$ = %0.0f' %xC[1]     + '   x$_e$ = %0.0f' %xC[0])
fig3.tight_layout() 


#%%  Fig 4: r =+10  t vs x
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (7,3)
fig4, ax = plt.subplots(1,3)

r = 10
xC = np.zeros(2); xC[0] = 0; xC[1] = r
C = 0
x0 = 12
tMax = 1
t = linspace(0,tMax,999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
ax[C].set_xlabel('t',color = 'black')
ax[C].set_ylabel('x',color= 'black')
ax[C].set_title(r'x(0) = %2.2f' %x0, fontsize = 12)
ax[C].grid()
ax[C].plot(t,xS,'b',lw = 2)    

C = 1
x0 = 5
tMax = 2.1
t = linspace(0,tMax,999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
ax[C].set_xlabel('t',color = 'black')
ax[C].set_ylabel('x',color= 'black')
ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
ax[C].grid()
ax[C].plot(t,xS,'b',lw = 2) 

C = 2
x0 = -0.1
tMax = 0.45
t = linspace(0,tMax,999)
sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 
ax[C].set_xlabel('t',color = 'black')
ax[C].set_ylabel('x',color= 'black')
ax[C].set_title('x(0) = %2.2f' %x0, fontsize = 12)
ax[C].grid()
ax[C].plot(t,xS,'r',lw = 2) 

fig4.suptitle('r = %0.0f' %r +'   x$_e$ = %0.0f' %xC[1]     + '   x$_e$ = %0.0f' %xC[0])
fig4.tight_layout() 


#%% r vs xE
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,4)
fig5, ax = plt.subplots(1)

ax.set_xlabel('r',color = 'black',fontsize = 12)
ax.set_ylabel('$x_C$',color= 'black',fontsize = 14)

ax.grid()

ax.set_xlim([-10,10]); ax.set_ylim([-10,10])

ax.plot([-10,0],[0,0],'b',lw = 2,label = 'Stable') 
ax.plot([0,10],[0,0],'r',lw = 2,label = 'Ustable') 
ax.plot([-10,0],[-10, 0],'r',lw = 2) 
ax.plot([0,10],[0, 10],'b',lw = 2) 

ax.legend() 
ax.axis('equal')
fig5.tight_layout() 

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

#%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')



'''




