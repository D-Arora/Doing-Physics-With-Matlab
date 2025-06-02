# -*- coding: utf-8 -*-
"""
S003.py               28 May 2025
 
 RAY OPTICS
 Matrix methods in paraxial optics
 Thick and thin lens matrices
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/S002.pdf

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
    A = 1; B = 0; C = 2/R; D = 1
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM4(R,n):        # M4 Refraction at a spherical interface 
    A = 1; B = 0; C = (n[0]-n[1])/(n[1]*R); D = n[0]/n[1]
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M

def TM5(R,n):        # M5 Thin len matrix 
    f = (n[1]-n[0])/n[0]*(1/R[0] - 1/R[1]); f = -1/f
                          
    A = 1; B = 0; C = -1/f; D = n[0]/n[1]
    M = zeros([2,2])
    M[0,0] = A; M[0,1] = B; M[1,0] = C; M[1,1] = D
    return M


#%% 1   Translation: translation matrix M
# Start point (0)
y0 = 5          # <<< y coordinate
a0Deg = 18      # <<< slope angle [deg]
L = 2           # <<< Translation distance
a0 = a0Deg*pi/180   # slope angle [rad]
# End point  (1)
V0 = vec(y0,a0)
M = TM1(L)
V1 = M@V0
y1 = V1[0,0]       # y coordinate
a1 = V1[1,0]       # slope angle [rad]
a1Deg = a1*180/pi  # slope angle [deg]
print('  ')
print('1   Translation')
print('     y0 = %0.3f' % y0 + '   a0 =  %0.3f  deg' % a0Deg)
print('     y1 = %0.3f' % y1 + '   a1 =  %0.3f  deg' % a1Deg)


#%% 2   Refraction at a plane interface
n = zeros(2)
y0 = 5          # <<< y coordinate
a0Deg = 18      # <<< slope angle [deg]

n[0] = 1       # <<< refractive index medium 0
n[1] = 1.5     # <<< refractive index medium 1

V0 = vec(y0,a0)
M = TM2(n)
V1 = M@V0
y1 = V1[0,0]       # y coordinate
a1 = V1[1,0]       # slope angle [rad]
a1Deg = a1*180/pi  # slope angle [deg]
print('  ')
print('2   Refraction at a plane interface')
print('     y0 = %0.3f' % y0 + '   a0 =  %0.3f  deg' % a0Deg)
print('     y1 = %0.3f' % y1 + '   a1 =  %0.3f  deg' % a1Deg)

 
#%% 3   Reflection from a spherical mirror
n = zeros(2)
y0 = 5          # <<< y coordinate
a0Deg = 18      # <<< slope angle [deg]
# radius of curvature: concave R < 0  /  convex R > 0  >>>
R = 10          

V0 = vec(y0,a0)
M = TM3(R)
V1 = M@V0
y1 = V1[0,0]       # y coordinate
a1 = V1[1,0]       # slope angle [rad]
a1Deg = a1*180/pi  # slope angle [deg]
print('  ')
print('3   Reflection from a spherical mirror')
print('     y0 = %0.3f' % y0 + '   a0 =  %0.3f  deg' % a0Deg)
print('     y1 = %0.3f' % y1 + '   a1 =  %0.3f  deg' % a1Deg)


#%% 4   Refraction at a spherical interface
n = zeros(2)
y0 = 5          # <<< y coordinate
a0Deg = 18      # <<< slope angle [deg]

# radius of curvature: concave R < 0  /  convex R > 0  >>>
R = 10 
n[0] = 1       # <<< refractive index medium 0
n[1] = 1.5     # <<< refractive index medium 1

V0 = vec(y0,a0)
M = TM4(R,n)
V1 = M@V0
y1 = V1[0,0]       # y coordinate
a1 = V1[1,0]       # slope angle [rad]
a1Deg = a1*180/pi  # slope angle [deg]
print('  ')
print('4  Refraction at a spherical interface')
print('     y0 = %0.3f' % y0 + '   a0 =  %0.3f  deg' % a0Deg)
print('     y1 = %0.3f' % y1 + '   a1 =  %0.3f  deg' % a1Deg)


#%% 5   Thin lens matrix 
R = zeros(2); n = zeros(2)
y0 = 5          # <<< y coordinate
a0Deg = 18      # <<< slope angle [deg]

# radius of curvature: >>>
R[0] = 10      # <<< radius of curvature 0
R[1] = 8     # <<< radius of curvature 1
n[0] = 1       # <<< refractive index medium 0
n[1] = 1.5     # <<< refractive index lens

V0 = vec(y0,a0)
M = TM5(R,n)
V1 = M@V0
y1 = V1[0,0]       # y coordinate
a1 = V1[1,0]       # slope angle [rad]
a1Deg = a1*180/pi  # slope angle [deg]
print('  ')
print('5  Thin lens matrix')
print('     y0 = %0.3f' % y0 + '   a0 =  %0.3f  deg' % a0Deg)
print('     y1 = %0.3f' % y1 + '   a1 =  %0.3f  deg' % a1Deg)


#%%  Example 1: propagation into long cylinder
L = 16     # object position
R = 4      # radius of curvature of spherical end of rod
n = zeros(2); n[0] = 1; n[1] = 1.5

T0 = TM1(L)       # Translation matrix: object
N  = TM4(R,n)     # Refraction matrix at spherical interface

y0 = 1
a0Deg = 4
a0 = a0Deg*pi/180
dy = L*a0
y1 = y0+dy
print(y1)

V0 = vec(y0,a0)
M = TM1(L)
V1 = M@V0
y1 = V1[0,0]       # y coordinate
a1 = V1[1,0]       # slope angle [rad]
a1Deg = a1*180/pi  # slope angle [deg]
print('  ')
print('1   Translation')
print('     y0 = %0.3f' % y0 + '   a0 =  %0.3f  deg' % a0Deg)
print('     y1 = %0.3f' % y1 + '   a1 =  %0.3f  deg' % a1Deg)

V0 = vec(y1,a1)
M = TM4(R,n)
V2 = M@V0
y2 = V2[0,0]       # y coordinate
a2 = V2[1,0]       # slope angle [rad]
a2Deg = a2*180/pi  # slope angle [deg]
print('  ')
print('4  Refraction at a spherical interface')
print('     y0 = %0.3f' % y1 + '   a0 =  %0.3f  deg' % a1Deg)
print('     y1 = %0.3f' % y2 + '   a1 =  %0.3f  deg' % a2Deg)


L = 24
V0 = vec(y2,a2)
M = TM1(L)
V3 = M@V0
y3 = V3[0,0]       # y coordinate
a3 = V3[1,0]       # slope angle [rad]
a3Deg = a3*180/pi  # slope angle [deg]
print('  ')
print('1   Translation')
print('     y0 = %0.3f' % y2 + '   a0 =  %0.3f  deg' % a2Deg)
print('     y1 = %0.3f' % y3 + '   a1 =  %0.3f  deg' % a3Deg)


#%%
L = 16
y0 = 1
a0Deg = 4
a0 = a0Deg*pi/180
dy = L*a0
y1 = y0+dy
V0 = vec(y0,a0)
M1 = TM1(L)
M2 = TM4(R,n)
x = 24
M3 = TM1(x)
M = M3@M2@M1

V = M@V0

print(V)
#x = 24           # unknown image position 
#print(T0)
#print(N)



#%% GRAPHICS
# # fig1, ax = plt.subplots(nrows=1, ncols=1)
# # plt.rcParams['font.size'] = 12
# # plt.rcParams["figure.figsize"] = (4,3)


# # ax.plot(xC,yC,'k',lw = 4)                          # circle

# # X = [0,xP]; Y = [0,yP]; ax.plot(X,Y,'k',lw = 1)    # radius vector 0P
# # X = [xA,xP]; Y = [yA,yP]; ax.plot(X,Y,'b',lw = 2)  # incident ray PA
# # X = [xB,xP]; Y = [yB,yP]; ax.plot(X,Y,'r',lw = 2)  # reflected ray BP

# # ax.plot(0,0,'ko',ms = 4)                           # Origin
# # ax.plot(abs(R/2),0,'ko',ms = 4)                    # Focal point

# # ax.set_xlim([-1,Ly+2]); ax.set_ylim([-6,6])
# # ax.set_xticks(arange(-2,Ly+0.5,2)); ax.set_yticks(arange(-6,6+0.5,2))

# # ax.set_xlabel('x', fontsize=12)
# # ax.set_ylabel('y', fontsize=12)
# # #ax.set_title('paraxial approx for slope  AE = %0.1f deg' %AE,fontsize = 10)
# # #ax.legend(fontsize = 10)
# # ax.grid()
# # ax.set_aspect('equal', 'box')

# # fig1.tight_layout()

# # #%% CONSOLE OUPUT
# # print('  ')
# # print('A(xA,yA) = (%0.2f' % xA + ', %0.2f)' %yA)
# # print('P(xP,yP) = (%0.2f, ' % xP + '%0.2f )' %yP)
# # print('B(xB,yB) = (%0.2f, ' % xB + '%0.2f )' %yB)
# # print('incidence slope angle   alpha0,  aI = %0.3f deg' %aIdeg)
# # print('reflection slope angle  alpha1, aR = %0.3f deg' %aRdeg)
# # print('angle of incidence   thetaI = %0.3f deg' %thetaIdeg)
# # print('angle of reflection  thetaR = %0.3f deg' %thetaRdeg)
# # print('radius vector phi = %0.3f deg' %phideg)      
# # print('focal length f = %0.3f ' % f)      

      
# # #%%
# # fig1.savefig('a1.png')      
      


