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


#%% INPUTS
# Refractive indices
n = 1; nL = 1.6
# Radii of curvature / thickness t
# CONCAVE R <0;  CONVEX R > 0
R1 = 45; R2 = -30; t = 5

#%% LENS MATRICES
M1 = zeros([2,2]); M2 = zeros([2,2]); M3 = zeros([2,2])

M00 = 1; M01 = 0; M10 = (n-nL)/(nL*R1); M11=n/nL 
M1[0,0] = M00; M1[1,0] = M01; M1[1,0] = M10; M1[1,1] = M11

M00 = 1; M01 = t; M10 = 0; M11=1
M2[0,0] = M00; M2[0,1] = M01; M2[1,0] = M10; M2[1,1] = M11

M00 = 1; M01 = 0; M10 = (nL-n)/(n*R2); M11=nL/n
M3[0,0] = M00; M3[0,1] = M01; M3[1,0] = M10; M3[1,1] = M11

M = M3@M2@M1

# focal length

f1= (nL-n)/n; f2 = (1/R1 - 1/R2); f3 = f1*f2
f = 1/f3


# N = 999
# R = -10              # radius of circle: *** concave surface R < 0  ***
# C1 = -45; C2 = 45  # circle limit angles   [deg]
# C = linspace(C1,C2,N); CR = C*pi/180
# xC = abs(R)*cos(CR); yC = abs(R)*sin(CR)
# Ly = 10;   # Y limit for grid

# # Points A, B and P
# xA = 0; yA = -3; xB = xA            # A and B
# yP = -1; xP = sqrt(R**2 - yP**2)    # P
# aI = (yP-yA)/(xP-xA)                 # alpha: slope incident ray AP

# phi = arctan(yP/xP)
#                              # angle alpha
# # Reflection matrix
# M = zeros([2,2]); V0 = zeros([2,1])
# M[0,0] = 1; M[1,0] = 2/R; M[1,1] = 1
# V0[0,0] = yP; V0[1,0]= aI
# V = M@V0
# aR = V[1,0]
# dx = xP - xA; dy = dx*aR
# yB = yP+dy 

# thetaI = aI-phi; thetaR = aR + phi
# q = 180/pi
# aIdeg = aI*q; aRdeg = aR*q
# thetaIdeg = thetaI*q; thetaRdeg = thetaR*q; phideg = phi*q

# # focal length
# f = abs(R) - yB/aR


# #%% GRAPHICS
# fig1, ax = plt.subplots(nrows=1, ncols=1)
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (4,3)


# ax.plot(xC,yC,'k',lw = 4)                          # circle

# X = [0,xP]; Y = [0,yP]; ax.plot(X,Y,'k',lw = 1)    # radius vector 0P
# X = [xA,xP]; Y = [yA,yP]; ax.plot(X,Y,'b',lw = 2)  # incident ray PA
# X = [xB,xP]; Y = [yB,yP]; ax.plot(X,Y,'r',lw = 2)  # reflected ray BP

# ax.plot(0,0,'ko',ms = 4)                           # Origin
# ax.plot(abs(R/2),0,'ko',ms = 4)                    # Focal point

# ax.set_xlim([-1,Ly+2]); ax.set_ylim([-6,6])
# ax.set_xticks(arange(-2,Ly+0.5,2)); ax.set_yticks(arange(-6,6+0.5,2))

# ax.set_xlabel('x', fontsize=12)
# ax.set_ylabel('y', fontsize=12)
# #ax.set_title('paraxial approx for slope  AE = %0.1f deg' %AE,fontsize = 10)
# #ax.legend(fontsize = 10)
# ax.grid()
# ax.set_aspect('equal', 'box')

# fig1.tight_layout()

# #%% CONSOLE OUPUT
# print('  ')
# print('A(xA,yA) = (%0.2f' % xA + ', %0.2f)' %yA)
# print('P(xP,yP) = (%0.2f, ' % xP + '%0.2f )' %yP)
# print('B(xB,yB) = (%0.2f, ' % xB + '%0.2f )' %yB)
# print('incidence slope angle   alpha0,  aI = %0.3f deg' %aIdeg)
# print('reflection slope angle  alpha1, aR = %0.3f deg' %aRdeg)
# print('angle of incidence   thetaI = %0.3f deg' %thetaIdeg)
# print('angle of reflection  thetaR = %0.3f deg' %thetaRdeg)
# print('radius vector phi = %0.3f deg' %phideg)      
# print('focal length f = %0.3f ' % f)      

      
# #%%
# fig1.savefig('a1.png')      
      


