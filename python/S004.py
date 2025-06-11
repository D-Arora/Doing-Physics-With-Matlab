# -*- coding: utf-8 -*-
"""
S004.py               06 June 2025
 
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
from numpy import arange, exp, linspace, zeros, amax, sqrt
from numpy import pi, sin, cos, tan, arctan
import matplotlib.pyplot as plt
import time

tStart = time.time()

plt.close('all')

#%%  INITIALIZE VARIABLES
#   a elevation (slope angle) [rad), aD slope angle [deg], y y-coordinate 
a = zeros(5); aD = zeros(5); y = zeros(5)
RD = 180/pi; DR = pi/180     # rad --> deg   deg --> rad
p = 0; q = 1                 # p input  /  q output

#%%  Inputs >>>
# Refractive indices
n0 = 1; n1 = 1.5 
# Convex front surface of lens: radius of curvature   R0 > 0 
R0 = 4   
# focal length: used to determine radius of concave rear of lens R1 < 0
f = 5  
# y coordinate point 0             
y[p] = 1                  
# Elevation (slope angle) 0 --> 1  [deg]
aD[p] = 12             # <<< slope angle [deg]
# Translation object plane to lens  0 --> 1: L01
L01 = 3.6
# Translation lens to image plane 3 --> 4: L34           
L34 = 10          
if L01 < f: L34 = f*L01/(L01-f)
# Concave rear surface of lens: radius of curvature   R1 > 0 
R1 = 1/R0 - n0/(f*(n1-n0))
R1 = 1/R1


#%% FUNCTIONS

def TM1(L):        # M1 Translation matrix 
    A = 1; B = L; C = 0; D = 1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM4(R,n0,n1):        # M4 Refraction at a spherical interface 
    A = 1; B = 0; C = (n0-n1)/(n1*R); D = n0/n1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M


#%%
def Tsetup(y0,aD0,L):
    a0 = DR*aD0
    M = TM1(L)
    y1,a1 = M@[y0,a0]
    aD1 = RD*a1
    return y1, aD1, M

def Rsetup(y0,aD0,R,n0,n1):
    a0 = DR*aD0
    M = TM4(R,n0,n1)
    y1,a1 = M@[y0,a0]
    aD1 = RD*a1
    return y1, aD1, M 

   
#%% Do transformations
# Translation 0 --> 1
y[q], aD[q], M1 = Tsetup(y[p],aD[p],L01)

# Refraction 1 --> 2
p = 1; q = 2
y[q], aD[q], M2 = Rsetup(y[p],aD[p],R0,n0,n1)

# Refraction 2 --> 3
p = 2; q = 3
y[q], aD[q], M3 = Rsetup(y[p],aD[p],R1,n1,n0)

# Translation 3 --> 4
p = 3; q = 4
y[q], aD[q], M4 = Tsetup(y[p],aD[p],L34)

a = aD*DR     # deg to rad

# System mmatrix
M = M4@M3@M2@M1
y1,a1 = M@[y[0],a[0]]
a1D = RD*a1
A = M[0,0]; B = M[0,1]; C = M[1,0]; D = M[1,1]


#%% Line from (x0,y0) through Origin (0,) / intersection point = image point
if L01 > f:
   xC = linspace(-L01,L34,999)
else:  
   xC = linspace(L34,20,999)
mC = -y[0]/L01
yC = mC*xC

# line 3 --> 4
bL = y[3]            # intercept
mL = (a[3])          # slope
x4 = bL/(mC - mL)    # intersection point x4
y4 = mL*x4+bL        # intersection point y4

xO = linspace(-20,20,999)
yO = a[3]*xO + y[3] 

#xxx

#%% GRAPHICS
plt.rcParams["figure.figsize"] = (6,3)
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1)

# Ray through centre of lens and output point
ax.plot(xC,yC,'m',lw = 1)  
ax.plot(xC[-1],yC[-1],'mo',ms = 6)

# Focal points 1 and 2
ax.plot(-f,0,'ko',ms = 6)    
ax.plot(f,0,'ko',ms = 6)      

# Optical axis and lens plane
X = [-20,26]; Y = [0,0]; ax.plot(X,Y,'k',lw = 1) 
X = [0,0]; Y = [-4,4]; ax.plot(X,Y,'k',lw = 3) 
# Image plane
X = [x4,x4]; Y = [-4,4]; ax.plot(X,Y,'k',lw = 1) 


# Input ray 0--> 1, input point (object), input plane
X = [-L01,0]; Y = [y[0],y[1]]; ax.plot(X,Y,'b',lw = 2) 
ax.plot(-L01,y[0],'bo',ms = 6)      
X = [-L01,-L01]; Y = [-4,4]; ax.plot(X,Y,'b',lw = 0.5) 

# Output ray 3 --> 4, output point (image), output plane
if L01 > f:
   X = [0,L34]; Y = [y[3],y[4]]; ax.plot(X,Y,'r',lw = 2)  
   ax.plot(L34,y[4],'ro',ms = 6)    
if L01 < f:     
   X = xO; Y = yO; ax.plot(X,Y,'r',lw = 2)  
   ax.plot(xC,yC,'m',lw = 1) 
   ax.plot(x4,y4,'ro',ms = 6) 
X = [L34,L34]; Y = [-4,4]; ax.plot(X,Y,'r',lw = 0.5) 
   


ax.set_xlim([-20,20]); ax.set_ylim([-4,4])
ax.set_title(r'f = %0.2f' %f +
     '  x$_0$ = %0.2f' %L01 + '  y$_0$ = %0.2f' %y[0] + 
     '  x$_4$ = %0.2f' %L34 + '  x$_I$ = %0.2f' %x4, fontsize = 10)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.grid()
fig1.tight_layout()

#xxx

#%% OUTPUT TO CONSOLE 
# print(' ')
# print('Propagation of a ray through a thin lens oprical system')
print('INPUTS')
print('   refractive index of air, n0 = %0.2f' % n0)
print('   refractive index of lens n1 = %0.2f' % n1)
print('   focal length, f = %0.2f' % f)
print('   radius convex front spherical surface of lens, R0 = %0.2f' % R0)
print('   Object plane,  L01 = %0.2f' % L01)
print('   Object height, y0 = %0.2f' % y[0])
print('   Incident ray elevation, a01 = %0.2f deg' % aD[0])
print('   Translation distance 3 --> 4,  L34 = %0.2f' % L34)
print('OUTPUTS  ')
print('   radius concave rear spherical surface of lens, R1 = %0.2f' % R1)
print('   refracted ray elevation, a34 = %0.2f' % aD[4])
print('   image plane,   x4 = %0.2f' % x4)
print('   image height,  y4 = %0.2f' % y4)
q = y[4]/y[0]; print('   lateral magnification, m = %0.2f' % q)
if a[0]!= 0:
   q = a[4]/a[0]
   print('   angular magnification, mA = %0.2f' % q)
print('   M --> A = %0.2f' % A + '  B = %0.2f' % B +
      '  C = %0.2f' % C + '  B = %0.2f' %D)
tExe = time.time() - tStart
print('\nExecution time %0.0f s' %tExe)

#%% Save figures
# fig1.savefig('a1.png')   
     
      


