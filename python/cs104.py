# -*- coding: utf-8 -*-
# cs104.py          Feb 2024

# NONLINEAR [1D] DYNAMICAL SYSTEMS
# FIXED POINTS, STABILITY ANALYSIS, BIFURCATIONS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_101.pdf

# https://math.libretexts.org/Bookshelves/Scientific_Computing_Simulations_and_Modeling/Scientific_Computing_(Chasnov)/II%3A_Dynamical_Systems_and_Chaos/12%3A_Concepts_and_Tools

# LIBRARIES  ================================================================
import numpy as np
from numpy import pi, sin, cos, linspace
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

tStart = time.time()


#%%   FUNCTIONS
def funct(r,x):
    xDot = r*x + x**3
    return xDot

#%%  SOLVE ODE for x
def lorenz(t, state):    
    xS = state
    dx = r*xS + xS**3
    return dx  


#%%  Fig 1:  Bifurcation diagrams
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7)
fig1, axes = plt.subplots(nrows=2, ncols=2)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.13,\
                    right = 0.95, hspace = 0.5,wspace=0.40)

# r = 0  ----------------------------------------------------------------------    
r = 0
N = 599; x1 = -6; x2 = 6; x = linspace(x1,x2,N)
y = funct(r,x)

R = 0; C = 0   
axes[R,C].set_ylabel('xDot',color= 'black',fontsize = 12)
axes[R,C].set_xlabel('x',color = 'black',fontsize = 12)
axes[R,C].set_title('r = %2.1f' % r, fontsize = 14)
axes[R,C].set_xlim([-6, 6])
axes[R,C].set_ylim([-200, 200])
axes[R,C].set_xticks(np.arange(-6,7,2))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
   
axes[R,C].plot(x, y, 'blue')
axes[R,C].plot(0, 0, 'ko',ms=8)

# r > 0  ----------------------------------------------------------------------
r = 16
y = funct(r,x)

R = 0; C = 1    
axes[R,C].set_ylabel('xDot',color= 'black',fontsize = 12)
axes[R,C].set_xlabel('x',color = 'black',fontsize = 12)
axes[R,C].set_title('r = %2.1f' % r, fontsize = 14)
axes[R,C].set_xlim([-6, 6])
axes[R,C].set_ylim([-200, 200])
axes[R,C].set_xticks(np.arange(-6,7,2))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
   
axes[R,C].plot(x, y, 'blue')
axes[R,C].plot(0, 0, 'ko',ms=8)

# r < 0  ----------------------------------------------------------------------
r = -16
x = linspace(-6,6,N)
y = funct(r,x)

R = 1; C = 0    
axes[R,C].set_ylabel('xDot',color= 'black',fontsize = 12)
axes[R,C].set_xlabel('x',color = 'black',fontsize = 12)
axes[R,C].set_title('r = %2.1f' % r, fontsize = 14)
axes[R,C].set_xlim([-6, 6])
axes[R,C].set_ylim([-100, 100])
axes[R,C].set_xticks(np.arange(-6,7,2))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()

axes[R,C].plot(x, y, 'blue')
axes[R,C].plot(0, 0, 'bo',ms=8)
axes[R,C].plot((-r)**0.5, 0, 'ro',ms=8)
axes[R,C].plot(-(-r)**0.5, 0, 'ro',ms=8)

# Bifurcation diagram
R = 1; C = 1    
axes[R,C].set_ylabel('$x_e$',color= 'black',fontsize = 14)
axes[R,C].set_xlabel('r',color = 'black',fontsize = 14)
axes[R,C].set_title('  ' , fontsize = 14) 
axes[R,C].set_xlim([-16, 16])
axes[R,C].set_xticks(np.arange(-16,17,4))

axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
 
r = linspace(-16,0,N)
xe = np.zeros(N)
axes[R,C].plot(r, xe, 'b',lw = 2)

r = linspace(0,16,N)
xe = np.zeros(N)
axes[R,C].plot(r, xe, 'r',lw = 2)

r = linspace(-16,0,N)
xe = (-r)**0.5
axes[R,C].plot(r, xe, 'r',lw = 2)
axes[R,C].plot(r, -xe, 'r',lw = 2)

fig1.savefig('a1.png')


#%% TIME EVOLUTION

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6,7)
fig1, axes = plt.subplots(nrows=3, ncols=1)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.13,\
                    right = 0.95, hspace = 0.5,wspace=0.40)

# r = 0  ----------------------------------------------------------------------    
r = 0
x0 = -1
N = 9999; t1 = 0; t2 = 0.499
t = linspace(t1,t2,N)

sol = odeint(lorenz, x0, t, tfirst=True)
x1 = sol[:,0] 

x0 = 1
sol = odeint(lorenz, x0, t, tfirst=True)
x2 = sol[:,0] 

R = 0  
axes[R].set_ylabel('xDot',color= 'black',fontsize = 12)
axes[R].set_xlabel('x',color = 'black',fontsize = 12)
axes[R].set_title('r = %2.1f' % r, fontsize = 14)
#axes[R].set_xlim([0, 1000])
#axes[R].set_ylim([-200, 200])
#axes[R].set_xticks(np.arange(-6,7,2))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
   
axes[R].plot(t, x1, 'blue')
axes[R].plot(t, x2, 'r')

# r > 0  ----------------------------------------------------------------------    
r = 16
x0 = -1
N = 9999; t1 = 0; t2 = 0.088
t = linspace(t1,t2,N)

sol = odeint(lorenz, x0, t, tfirst=True)
x1 = sol[:,0] 

x0 = 1
sol = odeint(lorenz, x0, t, tfirst=True)
x2 = sol[:,0] 

R = 1  
axes[R].set_ylabel('xDot',color= 'black',fontsize = 12)
axes[R].set_xlabel('x',color = 'black',fontsize = 12)
axes[R].set_title('r = %2.1f' % r, fontsize = 14)
#axes[R].set_xlim([0, 0.5])
#axes[R].set_ylim([-200, 200])
#axes[R].set_xticks(np.arange(-6,7,2))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
   
axes[R].plot(t, x1, 'blue')
axes[R].plot(t, x2, 'r')

# r > 0  ----------------------------------------------------------------------    
r = -16

N = 9999

x0 = -3
t1 = 0; t2 = 0.5
t = linspace(t1,t2,N)
sol = odeint(lorenz, x0, t, tfirst=True)
x1 = sol[:,0] 

x0 = 3
sol = odeint(lorenz, x0, t, tfirst=True)
x2 = sol[:,0] 

x0 = 0.5
sol = odeint(lorenz, x0, t, tfirst=True)
x3 = sol[:,0] 

x0 = -0.5
sol = odeint(lorenz, x0, t, tfirst=True)
x4 = sol[:,0]


R = 2  
axes[R].set_ylabel('xDot',color= 'black',fontsize = 12)
axes[R].set_xlabel('x',color = 'black',fontsize = 12)
axes[R].set_title('r = %2.1f' % r, fontsize = 14)
axes[R].set_xlim([0, 0.5])
axes[R].set_ylim([-6, 6])
#axes[R].set_xticks(np.arange(-6,7,2))
axes[R].set_yticks(np.arange(-6,7,2))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
   
axes[R].plot(t, x1, 'blue')
axes[R].plot(t, x2, 'r')
# axes[R].plot(t, x3, 'm') 
# axes[R].plot(t, x4, 'm')
       
fig1.savefig('a2.png')





