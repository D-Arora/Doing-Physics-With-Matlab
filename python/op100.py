# -*- coding: utf-8 -*-
"""
op100.py               oct 2024
 
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
 Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/op100.pdf
"""

import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array, dot, cross,sqrt
from numpy import degrees
from numpy.linalg import eig, norm
from scipy.integrate import odeint, quad, dblquad, simps

import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import axes3d
import sympy as sym
import math
from math import atan2, asin, acos

tStart = time.time()

#%% Cell 1: degrees <---> radians
# Input angle in degrees --> output in radians
thetaD = 180
thetaR = np.radians(thetaD)
print('Angle = %0.3f deg' % thetaD + '  = %0.3f rad' % thetaR)

# Input angle in radianss --> output in degrees
thetaR = pi/4
thetaD = np.degrees(thetaR)
print('Angle = %0.3f deg' % thetaD + '  = %0.3f rad' % thetaR)

#%% Cell #2:  [3D] vector
V1 = array([3.0,5,-6]) 
V2 = array([5.0, 9, 2])
V3 = array([5.0, 9, 2])
V3.shape = (3,1)

# Addition, subtraction, magnitudes, azimuthal angle
V12 = V1 + V2
V21 = V1 - V2
V13 = V1 + V3
V31 = V1 - V3

V1mag = norm(V1)
V2mag = norm(V2)
V3mag = norm(V3)
V12mag = norm(V12)
V21mag = norm(V21)
V13mag = norm(V13)
V31mag = norm(V31)
phi2 = sqrt(V2[0]**2 + V2[1]**2)

# Return the arc tangent of y/x in radians: Azimuthal and polar angles
phi2R = atan2(V2[1],V2[0])
phi2Deg = np.degrees(phi2R)

#%% Cell #3: [3D] plot of a vector
x = V2[0]; y = V2[1]; z = V2[2]
plt.rcParams['font.size'] = 10
fig1 = plt.figure(figsize=(4,3))
ax = plt.axes(projection='3d')
X = [0,x]; Y = [0,y]; Z = [0,z]
ax.plot3D(X,Y,Z,'b',lw = 2)
ax.plot3D(x,y,z,'ro',ms = 6)
ax.plot3D(0,0,0,'ko')
ax.plot3D([0,x],[0,0],[0,0],'m')
ax.plot3D([0,0],[0,y],[0,0],'m')
ax.plot3D([0,0],[0,0],[0,z],'m')
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.set_zlabel('z', fontsize=12)
fig1.tight_layout()
ax.view_init(19,-103,0)
# fig1.savefig('a1.png') 

#%% Cell #4  Vector multiplication

M1 = V1*V2
M2 = V2*V1
M3 = V1*V3
M4 = V3*V1

# Dot product
mDot12 = np.dot(V1,V2)
mDot21 = np.dot(V2,V1)
mDot13 = np.dot(V1,V3)
# mDot31 = np.dot(V3,V1)    wrong dimensions

#%% Cell #5:  Array manipulations
A = array([1.0,0,1])
B = array([0.0,1,1])
C = array([2.0,4,5])

# Dot product
AdotB = dot(A,B)
BdotA = dot(B,A)
AdotC = dot(A,C)
CdotA = dot(C,A)
BdotC = dot(B,C)
CdotB = dot(C,B)
# Magnitudes and angles between vectors
Amag = norm(A)
Bmag = norm(B)
Cmag = norm(C)
ABangle = acos(AdotB/(Amag*Bmag))
ABdeg   = degrees(ABangle)
BAangle = acos(BdotA/(Amag*Bmag))
BAdeg   = degrees(BAangle)
ACangle = acos(AdotC/(Amag*Cmag))
ACdeg   = degrees(ACangle)
BCangle = acos(BdotC/(Bmag*Cmag))
BCdeg   = degrees(BCangle)
# Cross products and their magnitude
AB = cross(A,B)
ABmag = norm(AB)
BA = cross(B,A)
BAmag = norm(BA)
AC = cross(A,C)
ACmag = norm(AC)
CA = cross(C,A)
CAmag = norm(CA)
BC = cross(B,C)
BCmag = norm(BC)
CB = cross(C,B)
CBmag = norm(CB)
# Triple products
# AdotBC = dot(A,BC)
# BdotCA = dot(B,CA)
# CdotAB = dot(C,AB)
# AcrossBC = cross(A,cross(B,C))
#AcrossBC1 = dot(B,)

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
