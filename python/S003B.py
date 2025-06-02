# -*- coding: utf-8 -*-
"""
S003B.py               01 June 2025
 
 RAY OPTICS
 Matrix methods in paraxial optics
   TRANSFOMATION MATRICES: ABCD matrices
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/S003.pdf

Functions to generate the tranfomation matrices ABCD

Example 1   Propagation into long cylinde
   translation 0 --> 1
   refraction  1 --> 2
   translation 2 --> 3
"""

import numpy as np
from numpy import arange, exp, linspace, zeros, amax, sqrt
from numpy import pi, sin, cos, tan, arctan
import matplotlib.pyplot as plt
import time
from scipy import special
from scipy.signal import find_peaks
from scipy.special import jv
from mpl_toolkits.mplot3d import Axes3D
tStart = time.time()

#plt.close('all')


#%% FUNCTIONS
# def vec(u,v):       # column vector [u,v]
#     V = zeros([2,1])
#     V[0,0] = u; V[1,0] = v
#     return V

def TM1(L):        # M1 Translation matrix 
    A = 1; B = L; C = 0; D = 1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM2(n):        # M2 Refraction at a plane interface 
    A = 1; B = 0; C = 0; D = n[0]/n[1]
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM3(R):        # M3 Reflection from spherical mirror
    A = 1; B = 0; C = 2/R[0]; D = 1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM4(R,n):        # M4 Refraction at a spherical interface 
    A = 1; B = 0; C = (n[0]-n[1])/(n[1]*R[0]); D = n[0]/n[1]
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM5(R,n):        # M5 Thin len matrix 
    f = (n[1]-n[0])/n[0]*(1/R[0] - 1/R[1]); f = -1/f
                          
    A = 1; B = 0; C = -1/f; D = n[0]/n[1]
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

def Rsetup(y0,aD0,R):
    a0 = DR*aD0
    M = TM4(R,n)
    y1,a1 = M@[y0,a0]
    aD1 = RD*a1
    return y1, aD1, M 

   
#%%  INITIALIZE VARIABLES
#       n refractive indices, R radii
#       a slope angle [rad), aD slope angle [deg], y y-coordinate 
n = zeros(4); R = zeros(4); a = zeros(4)
aD = zeros(4); y = zeros(4)
RD = 180/pi; DR = pi/180     # rad --> deg   deg --> rad
p = 0; q = 1                 # p input  /  q output


#%%  Inputs >>>
n[0] = 1; n[1] = 1.5      # Refractive indices
R[0] = 4                  # Radius of curvature of spherical end of rod

y[p] = 1                  # <<< y coordinate

aD[p] = -5                # <<< slope angle [deg]

                          # Translations
L01 = 16                 # Object plane to cylinder
 
L23 = 20                  # <<< Cylinder to image plane


#%% Do transformations
# Translation 0 --> 1
y[q], aD[q], M1 = Tsetup(y[p],aD[p],L01)

# Refraction 1 --> 2
p = 1; q = 2
y[q], aD[q], M2 = Rsetup(y[p],aD[p],R)

# Translation 2 --> 3
p = 2; q = 3
y[q], aD[q], M3 = Tsetup(y[p],aD[p],L23)

a = aD*DR     # deg to rad

# System mmatrix
M = M3@M2@M1
y1,a1 = M@[y[0],a[0]]
a1D = RD*a1
A = M[0,0]; B = M[0,1]; C = M[1,0]; D = M[1,1]

# OUTPUT TO CONSOLE 
print(' ')
print('Example 1: propagation into long cylinder  ')
print('Object plane  L01 = %0.2f' % L01)
print('Image plane  x = L23 = %0.2f' % L23)
print('Input vector:  y0 = %0.3f' % y[0] + '   a0 =  %0.3f  deg' % aD[0])
print('0 --> 1   Translation')
print('   y1 = %0.3f' % y[1] + '   a1 =  %0.3f  deg' % aD[1])
print(M1)
print('')
print('1 --> 2   Refraction')
print('   y2 = %0.3f' % y[2] + '   a2 =  %0.3f  deg' % aD[2])
print(M2)
print('')
print('2--> 3   Translation')
print('   y3 = %0.3f' % y[3] + '   a3 =  %0.3f  deg' % aD[3])
print(M3)
print('')
print('System matrix')
print('Input vector:  y0 = %0.3f' % y[0] + '   a0 =  %0.3f  deg' % aD[0])
print('Output vector:  y1 = %0.3f' % y1 + '    a1 =  %0.3f  deg' % a1D)
print('A = %0.2f' %A + '   B = %0.2f' %B +
      '   C = %0.2f' %C + '   D = %0.2f' %D)

#xxx

#%% GRAPHICS
plt.rcParams["figure.figsize"] = (4,3)
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1)

X = [-L01,0]; Y = [y[0],y[1]]
ax.plot(X,Y,'b',lw = 2) 
ax.plot(-L01,1,'bo',ms = 6)  

X = [0,L23]; Y = [y[2],y[3]]
ax.plot(X,Y,'r',lw = 2)     
ax.plot(L23,y[3],'ro',ms = 6)  

X = [-20,26]; Y = [0,0]
ax.plot(X,Y,'k',lw = 1) 

ax.set_xlim([-20,26]); ax.set_ylim([-2,4])
ax.set_yticks(arange(-2,4.5,1))
#ax.set_title('y$_0$ = %0.2f' %y[0] + '   y$_1$ = %0.2f' %y[q] +
#              '  x = %0.2f' %L23,fontsize = 10)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.grid()
fig1.tight_layout()

# fig1.savefig('a1.png')   


      
      


