# -*- coding: utf-8 -*-
"""
S005.py               06 June 2025
 
 RAY OPTICS
   Matrix methods in paraxial optics
   TRANSFOMATION MATRICES: ABCD matrices
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/S003.pdf

Functions to generate the transformation matrices ABCD

THIN  LENS  --> optical path through the optical system
"""

import numpy as np
from numpy import arange, exp, linspace, zeros, amax, sqrt,array
from numpy import pi, sin, cos, tan, arctan
import matplotlib.pyplot as plt
import time

tStart = time.time()

plt.close('all')

#%%  INITIALIZE VARIABLES
#   optical axis x / height y, elevation (slope angle) a [rad), 
x = zeros(4); y = zeros([4,3]); a = zeros([4,3])


#%%  Inputs >>>
# focal length
f = 5  
# yI  coordinate point 0             
yI = 2                  
# Input elevation (slope angle)   [deg]
aIdeg = array([18,0,-18])  
# Translation object plane to lens  0 --> 1
sI = 3
# Translation lens to output plane2 --> 3        
sF = 20


#%%  SETUP
s1 = f*sI/(sI-f)
aI = aIdeg*pi/180
x[0]  = -sI
x[-1] = sF
y[0,:] = yI
a[0,:] = aI
xI = x[0]

if sI < f: sF = s1

#%% FUNCTIONS

def TM(L):        # M Translation matrix 
    A = 1; B = L; C = 0; D = 1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def RM(f):        # Refraction matrix thin lens 
    A = 1; B = 0; C = -1/f; D = 1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M


#%%
def Tsetup(y0,a0,L):     # angle [rad]
    M = TM(L)
    y1,a1 = M@[y0,a0]
    return y1, a1, M

def Rsetup(y0,a0):# angle [rad]
    M = RM(f)
    y1,a1 = M@[y0,a0]
    return y1, a1, M 

   
#%% Do transformations
# Translation 
p = 0; q = 1
for c in range(3):
    Y = y[p,c]
    A = a[p,c]
    S = sI
    y[q,c], a[q,c], M1 = Tsetup(Y,A,S)
    
p = p+1; q = q+1
for c in range(3):
    Y = y[p,c]
    A = a[p,c]
    y[q,c], a[q,c], M2 = Rsetup(Y,A)

p = p+1; q = q+1
for c in range(3):
    Y = y[p,c]
    A = a[p,c]
    S = sF
    y[q,c], a[q,c], M1 = Tsetup(Y,A,S)

xF = x[-1]

#%% System mmatrix
M1 = TM(sI)
M2 = RM(f)
M3 = TM(sF)
M = M3@M2@M1

yF,aF = M@[y[0],a[0]]

aFdeg = aF*180/pi

A = M[0,0]; B = M[0,1]; C = M[1,0]; D = M[1,1]

# Intersection point of output rays (xQ,yQ) and magnification
xQ = zeros(2)
b0 = y[2,0]; m0 = a[2,0]
b1 = y[2,1]; m1 = a[2,1]
b2 = y[2,2]; m2 = a[2,2]
xQ[0] = (b2-b1)/(m1-m2); xQ[1] = (b2-b0)/(m0-m2)
yQ = m1*xQ+b1
mag = yQ/yI

# s0 < f   virtual image
XQ = xQ[0]; YQ = yQ[0]
x0 = linspace(XQ,20,999)
b = YQ - a[3]*XQ
y0 = a[3,0]*x0+b[0]
y1 = a[3,1]*x0+b[1]
y2 = a[3,2]*x0+b[2]

#%% GRAPHICS
plt.rcParams["figure.figsize"] = (6,3)
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1)

# Input ray 0--> 1, input point (object), input plane
for c in range(3):
    X = [x[0],x[1]]; Y = [y[0,0],y[1,c]]; ax.plot(X,Y,'b',lw = 1.5) 
    if sI > f:
       X = [x[2],x[3]]; Y = [y[2,c],y[-1,c]]; ax.plot(X,Y,'r',lw = 1.5) 
    if sI < f:
       X = x0; Y = y0; ax.plot(X,Y,'r',lw = 1.5)  
       Y = y1; ax.plot(X,Y,'r',lw = 1.5) 
       Y = y2; ax.plot(X,Y,'r',lw = 1.5)
       
# Focal points 1 and 2
ax.plot(-f,0,'ko',ms = 6)    
ax.plot(f,0,'ko',ms = 6)      

yLim = ax.get_ylim()
# Optical axis and lens plane
X = [-20,26]; Y = [0,0]; ax.plot(X,Y,'k',lw = 1) 
X = [0,0]; Y = yLim; ax.plot(X,Y,'k',lw = 3) 

# Output plane and focal plane
X = [s1,s1]; Y = yLim; ax.plot(X,Y,'m',lw = 1) 
X = [sF,sF]; Y = yLim; ax.plot(X,Y,'m',lw = 1) 

# Input output point
ax.plot(xI,yI,'bo',ms = 6)      

ax.set_xlim([-21,21])     
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.grid()
fig1.tight_layout()


#%% OUTPUT TO CONSOLE 
print(' ')
print('Propagation of a ray through a thin lens optical system')
print('focal length, f = %0.2f' % f)
print('Intersection point (image point)  xQ = %.2f' %xQ[0] + '   yQ = %0.2f' %yQ[0])
print('Magnification at image point mag = %0.2f ' %mag[0])
np.set_printoptions(precision=3)
print('position x')
print(x)
print('heights y')
print(y)
print('Elevations [deg]')
q = a*180/pi; print(q)
print('  ')
print('System matrix and Cardinal points') 
print('   M --> A = %0.2f' % A + '  B = %0.2f' % B +
      '  C = %0.2f' % C + '  D = %0.2f' %D)
F1 = xI+D/C; F2 = xF-A/C; print('   Focal distance F1 = %0.2f'%F1 +
            '   F2 = %0.2f' %F2) 
H1 = xI+(D-1)/C; H2 = xF+(1-A)/C;print('   Principle plane H1 = %0.2f'%H1 +
            '   H2 = %0.2f' %H2)
N1 = xI+(D-1)/C; N2 = xF+(1-A)/C;print('   Nodal plane N1 = %0.2f'%N1 +
            '   N2 = %0.2f' %N2)
 
tExe = time.time() - tStart
print('\nExecution time %0.0f s' %tExe)

#%% Save figures
# fig1.savefig('a1.png')   
     
      


