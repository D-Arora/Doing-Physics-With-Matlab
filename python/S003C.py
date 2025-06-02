# -*- coding: utf-8 -*-
"""
S003C.py               01 June 2025
 
 RAY OPTICS
 Matrix methods in paraxial optics
   TRANSFOMATION MATRICES: ABCD matrices
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/S003.pdf

Functions to generate the tranfomation matrices ABCD

vec --> create a 2 row column vector [u,v]
TM1 --> Translation matrix 
TM2 --> Refraction at a plane interface  matrix 
TM3 --> Reflection from spherical mirror matrix
TM4 --> Refraction at a spherical interface  matrix 
TM5 --> Thin len matrix

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

plt.close('all')


#%% FUNCTIONS
def vec(u,v):       # column vector [u,v]
    V = zeros([2,1])
    V[0,0] = u; V[1,0] = v
    return V

def TM1(L):        # Translation matrix 
    A = 1; B = L; C = 0; D = 1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM2(n):        # Refraction at a plane interface 
    A = 1; B = 0; C = 0; D = n[0]/n[1]
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM3(R):        # Reflection from spherical mirror
    A = 1; B = 0; C = 2/R[0]; D = 1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM4(R,n):        # Refraction at a spherical interface 
    A = 1; B = 0; C = (n[0]-n[1])/(n[1]*R[0]); D = n[0]/n[1]
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM5(R,n):        # Thin len matrix 
    f = (n[1]-n[0])/n[0]*(1/R[0] - 1/R[1]); f = -1/f
                          
    A = 1; B = 0; C = -1/f; D = n[0]/n[1]
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M


#%%  INITIALIZE VARIABLES
#       n refractive indices, R radii
#       a slope angle [rad), aD slope angle [deg], y y-coordinate 
n = zeros(4); R = zeros(4); a = zeros(4)
aD = zeros(4); y = zeros(4)
RD = 180/pi; DR = pi/180     # rad --> deg   deg --> rad

#%% 1   Translation: translation matrix M
p = 0   # start point
q = 1   # end point
# Start point (0)
L = 2          # <<< translation distance 
y[p] = 5          # <<< y coordinate
aD[p] = 18        # <<< slope angle [deg]
a[p] = aD[p]*DR      # slope angle [rad]
u = y[p]; v = a[p] 
V0 = vec(u,v)
M = TM1(L)
V1 = M@V0
y[q] = V1[0,0]       # y coordinate
a[q] = V1[1,0]       # slope angle [rad]
aD[q] = a[q]*RD         # slope angle [deg]
print('  ')
print('1   Translation')
print('     y0 = %0.3f' % y[p] + '   a0 =  %0.3f  deg' % aD[q])
print('     y1 = %0.3f' % y[q] + '   a1 =  %0.3f  deg' % aD[q])


#%% 2   Refraction at a plane interface
p = 0   # start point
q = 1   # end point
n[p] = 1       # <<< refractive index medium 0
n[q] = 1.5     # <<< refractive index medium 1

# Start point (0)
y[p] = 5          # <<< y coordinate
aD[p] = 18        # <<< slope angle [deg]
a[p] = aD[p]*DR      # slope angle [rad]
u = y[p]; v = a[p] 
V0 = vec(u,v)
M = TM2(n)
V1 = M@V0
y[q] = V1[0,0]       # y coordinate
a[q] = V1[1,0]       # slope angle [rad]
aD[q] = a[q]*RD      # slope angle [deg]

print('  ')
print('2   Refraction at a plane interface')
print('     y0 = %0.3f' % y[p] + '   a0 =  %0.3f  deg' % aD[p])
print('     y1 = %0.3f' % y[q] + '   a1 =  %0.3f  deg' % aD[q])

 
#%% 3   Reflection from a spherical mirror
p = 0   # start point
q = 1   # end point
R[p] = 10     # <<< refractive index medium 1

# Start point (0)
y[p] = 5             # <<< y coordinate
aD[p] = 18           # <<< slope angle [deg]
a[p] = aD[p]*DR      # slope angle [rad]
u = y[p]; v = a[p] 
V0 = vec(u,v)
M = TM3(R)
V1 = M@V0
y[q] = V1[0,0]       # y coordinate
a[q] = V1[1,0]       # slope angle [rad]
aD[q] = a[q]*RD      # slope angle [deg]

print('  ')
print('3   Reflection from a spherical mirror')
print('     y0 = %0.3f' % y[p] + '   a0 =  %0.3f  deg' % aD[p])
print('     y1 = %0.3f' % y[q] + '   a1 =  %0.3f  deg' % aD[q])


#%% 4   Refraction at a spherical interface
p = 0   # start point
q = 1   # end point
n[p] = 1       # <<< refractive index medium 0
n[q] = 1.5     # <<< refractive index medium 1
R[p] = 10      # <<<  radius
# Start point (0)
y[p] = 5             # <<< y coordinate
aD[p] = 18           # <<< slope angle [deg]
a[p] = aD[p]*DR      # slope angle [rad]
u = y[p]; v = a[p] 
V0 = vec(u,v)
M = TM4(R,n)
V1 = M@V0
y[q] = V1[0,0]       # y coordinate
a[q] = V1[1,0]       # slope angle [rad]
aD[q] = a[q]*RD      # slope angle [deg]

print('  ')
print('4  Refraction at a spherical interface')
print('     y0 = %0.3f' % y[p] + '   a0 =  %0.3f  deg' % aD[p])
print('     y1 = %0.3f' % y[q] + '   a1 =  %0.3f  deg' % aD[q])


#%% 5   Thin lens matrix 
p = 0   # start point
q = 1   # end point

R[p] = 10      # <<< radius of curvature 0
R[q] = 8       # <<< radius of curvature 1
n[p] = 1       # <<< refractive index medium 0
n[q] = 1.5     # <<< refractive index lens

# Start point (0)
y[p] = 5             # <<< y coordinate
aD[p] = 18           # <<< slope angle [deg]
a[p] = aD[p]*DR      # slope angle [rad]
u = y[p]; v = a[p] 
V0 = vec(u,v)
M = TM5(R,n)
V1 = M@V0
y[q] = V1[0,0]       # y coordinate
a[q] = V1[1,0]       # slope angle [rad]
aD[q] = a[q]*RD      # slope angle [deg]

print('  ')
print('5  Thin lens matrix')
print('     y0 = %0.3f' % y[p] + '   a0 =  %0.3f  deg' % aD[p])
print('     y1 = %0.3f' % y[q] + '   a1 =  %0.3f  deg' % aD[q])




