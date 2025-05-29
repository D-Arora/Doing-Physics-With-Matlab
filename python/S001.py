# -*- coding: utf-8 -*-
"""
S001.py               24 May 2025
 
 RAY OPTICS
 Fermat’s Principle: calculus of variations - variational principle  
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/S001.pdf

"""

import numpy as np
from numpy import pi, arange, exp, sin, linspace, zeros, amax, sqrt, real
import matplotlib.pyplot as plt
import time
from scipy import special
from scipy.signal import find_peaks
from scipy.special import jv
from mpl_toolkits.mplot3d import Axes3D
tStart = time.time()

plt.close('all')

#%% REFLECTION FROM PLANE MIRROR
# INPUT PARAMETERS:    x < 0
x = zeros(3); y =zeros(3); n = zeros(3)

# Reflection inputs >>>
x[0] = -2;  y[0] = 2
x[1] = -2;  y[1] = -2

# Refraction inputs >>>
x[2] = 5;  y[2] = -5
n[2] = 1.5; n[0] = 1; n[1] = 1

Y1 = min(y); Y2 = max(y); N = 9999

# SETUP AND CALCULATIONS: REFLECTION
Y = linspace(Y1,Y2,N); dY = Y[2] - Y[1]    # position along Y-axis
L1 = sqrt((0-x[1])**2 + (Y-y[1])**2)       # length of incident ray
L0 = sqrt((0-x[0])**2 + (Y-y[0])**2)       # length of reflected ray
L01 = L0 + L1                                # total path length
m = np.gradient(L01,dY)             # slope m = dL/dY
Y01 = Y[L01==min(L01)][0]              # Y when L = min   point P
Y01m = Y[m**2 == min(m**2)][0]  # Y when dL/dY = 0
L01min = min(L01)                    # min total path length
# Angles of incidence (0) and reflection (1)   [deg]
theta = np.arctan((abs((y-Y01)/x))) * (180/pi)  

# GRAPHICS: REFLECTION
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax1 = ax.twinx()
ax.plot(Y,L01,'b',lw = 2)
ax1.plot(Y,m,'r',lw = 2)
ax1.plot([Y1,Y2],[0,0],'r',lw = 1)
ax.set_xlabel('Y [m]', fontsize=10)
ax.set_ylabel('L [m]', fontsize=10)
ax1.set_ylabel('dL/de', fontsize=1)
ax.set_title('P(y) = %0.2f m' %Y01)
ax.grid()
fig1.tight_layout()

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.plot([x[0],0],[y[0],Y01],'b')
ax.plot([x[1],0],[y[1],Y01],'r')
ax.plot([0,0],[Y1,Y2],'k',lw = 1)
ax.plot([min(x),0],[Y01,Y01],'k',lw = 1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal', 'box')
ax.set_title('$\Theta_0$ = %0.2f deg' %theta[0]
         + '  $\Theta_1$ = %0.2f deg' %theta[1],fontsize = 10 )
ax.set_aspect('equal', 'box')
ax.grid() 
fig2.tight_layout()



#%%
# SETUP AND CALCULATIONS: REFRACTION
Y = linspace(Y1,Y2,N); dY = Y[2] - Y[1]    # position along Y-axis
L2 = n[2]*sqrt((0-x[2])**2 + (Y-y[2])**2)       # length of incident ray
L0 = n[0]*sqrt((0-x[0])**2 + (Y-y[0])**2)       # length of reflected ray
L02 = L0 + L2                                # total path length
m = np.gradient(L02,dY)             # slope m = dL/dY
Y02 = Y[L02==min(L02)][0]          # Y when L = min
Y02m = Y[m**2 == min(m**2)][0]  # Y when dL/dY = 0
L02min = min(L02)                    # min total path length
# Angles of incidence (0) and reflection (1)   [deg]
theta[2] = np.arctan((abs((y[2]-Y02)/x[2]))) * (180/pi)  

# GRAPHICS  REFRACTION
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax1 = ax.twinx()

ax.plot(Y,L02,'b',lw = 2)
ax1.plot(Y,m,'r',lw = 2)
ax1.plot([Y1,Y2],[0,0],'r',lw = 1)

ax.set_xlabel('Y [m]', fontsize=10)
ax.set_ylabel('L [m]', fontsize=10)
ax1.set_ylabel('dL/de', fontsize=1)
ax.set_title('P(y) = %0.2f m' %Y02)
ax.grid()
fig3.tight_layout()

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,2.5)
fig4, ax = plt.subplots(nrows=1, ncols=1)

ax.plot([x[0],0],[y[0],Y02],'b')
d = abs(y[0]-Y02)
ax.plot([0,x[0]],[Y02,Y02-d],'b')
ax.plot([x[2],0],[y[2],Y02],'r')
ax.plot([0,0],[Y1,Y2],'k',lw = 1)
ax.plot([min(x),max(x)],[Y02,Y02],'k',lw = 1)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal', 'box')
ax.set_title('$\Theta_0$ = %0.2f deg' %theta[0]
         + '  $\Theta_2$ = %0.2f deg' %theta[2],fontsize = 10 )
ax.set_aspect('equal', 'box')
ax.grid() 
fig4.tight_layout()

# CONSOLE 
print('  ')
print('REFLECTION')
print('L_min = %0.3f  m' % L01min)
print('L = L_min   Y = %0.3f  m' % Y01)  
print('dL/de = 0   Y = %0.3f  m' % Y01m) 
print('theta0 = %0.3f  deg' % theta[0])   
print('theta1 = %0.3f  deg' % theta[1]) 

print('  ')
print('REFRACTION')
print('L_min = %0.3f  m' % L02min)
print('L = L_min   Y = %0.3f  m' % Y02)  
print('dL/de = 0   Y = %0.3f  m' % Y02m) 
print('theta0 = %0.3f  deg' % theta[0])   
print('theta2 = %0.3f  deg' % theta[2]) 
q = n[0]*sin(theta[0])
q = n[2]*sin(theta[2]) 
tExe = time.time() - tStart
print('\nExecution time %0.0f s' %tExe)

#%%

"""
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')

"""

