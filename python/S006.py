# -*- coding: utf-8 -*-
"""
S006.py               18 June 2025
 
 RAY OPTICS
   Matrix methods in paraxial optics
   TRANSFOMATION MATRICES: ABCD matrices
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/S003.pdf

Functions to generate the transformation matrices ABCD
3 optical paths through the optical system
Compound lens system: ltwo converging lens 


https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_III_-_Optics_and_Modern_Physics_(OpenStax)/02%3A_Geometric_Optics_and_Image_Formation/2.09%3A_Microscopes_and_Telescopes

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
D = 6
x = zeros(D); y = zeros([D,3]); a = zeros([D,3])


#%%  Inputs >>>
# focal lengths
fA = 5; fB = 8
# Translations
L01 = 8
L23 = 20
L45 = 20
# y coordinate point /  elevation (slope angle)   [deg]             
y[0,:] = 1.0                 
a0deg = array([18,0,-18])  

# X optical axis limits
X1 = -20; X2 = 40


#%%  lens A equation
a[0] = a0deg*pi/180
sA0 = L01
sA1 = fA*sA0/(sA0 - fA)
# x positions
x[0] = 0; x[1] = x[0] + L01; x[2] = x[1]
x[3] = x[2] + L23; x[4] = x[3]; x[5] = x[4] + L45   


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

def Tsetup(y0,a0,L):     # angle [rad]
    M = TM(L)
    y1,a1 = M@[y0,a0]
    return y1, a1, M

def Rsetup(f,y0,a0):# angle [rad]
    M = RM(f)
    y1,a1 = M@[y0,a0]
    return y1, a1, M 

   
#%% Do transformations
# Translation 0 --> 1  L01
p = 0; q = 1
for c in range(3):
    Y = y[p,c]
    A = a[p,c]
    S = L01
    y[q,c], a[q,c], M1 = Tsetup(Y,A,S)

# Refraction lens A    
p = p+1; q = q+1
for c in range(3):
    Y = y[p,c]
    A = a[p,c]
    y[q,c], a[q,c], M2 = Rsetup(fA,Y,A)

# Translation 2 --> 3  L23
p = p+1; q = q+1
for c in range(3):
    Y = y[p,c]
    A = a[p,c]
    S = L23
    y[q,c], a[q,c], M1 = Tsetup(Y,A,S)

# Refraction lens B    
p = p+1; q = q+1
for c in range(3):
    Y = y[p,c]
    A = a[p,c]
    y[q,c], a[q,c], M2 = Rsetup(fB,Y,A)

# Translation 4 --> 5  L45
p = p+1; q = q+1
for c in range(3):
    Y = y[p,c]
    A = a[p,c]
    S = L45
    y[q,c], a[q,c], M1 = Tsetup(Y,A,S)
    

#%% System mmatrix
TM1 = TM(L01)
RM1 = RM(fA)
TM2 = TM(L23)
RM2 = RM(fB)
TM3 = TM(L45)
M = TM3 @ RM2 @ TM2 @ RM1 @ TM1

yF,aF = M@[y[0],a[0]]

aFdeg = aF*180/pi

A = M[0,0]; B = M[0,1]; C = M[1,0]; D = M[1,1]


#%% Intersection point of output rays (xA,yA) (xB,yB) / magnifications mT mA
xB = zeros(2); bB = zeros(3); xA = zeros(2); bA = zeros(3)
for c in range(3):
    bB[c] = y[4,c] - a[4,c]*x[4]
xB[0] = (bB[2]-bB[1])/(a[4,1]-a[4,2])
xB[1] = (bB[2]-bB[0])/(a[4,0]-a[4,2])
yB = a[4,0]*xB+bB[0]
XB = linspace(-60,x[4],999)
YB0 = a[4,0]*XB + bB[0] 
YB1 = a[4,1]*XB + bB[1] 
YB2 = a[4,2]*XB + bB[2]

for c in range(3):
    bA[c] = y[2,c] - a[2,c]*x[2]
xA[0] = (bA[2]-bA[1])/(a[2,1]-a[2,2])
xA[1] = (bA[2]-bA[0])/(a[2,0]-a[2,2])
yA = a[2,0]*xA+bA[0]

# Transverse magnification
mT = yB[0]/y[0]
# Angular magnification
alpha0 = y[0][0] / x[2]
alpha1 = yB[0] / (x[4]-xB[0])
mA = alpha1 / alpha0


#%% GRAPHICS
plt.rcParams["figure.figsize"] = (6,3)
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1)

# Ray paths
for c in range(3):
    Y = y[:,c]
    ax.plot(x,Y,'k',lw = 1.5) 
           
# Focal points A and B
xA0 = x[1] - fA; ax.plot(xA0,0,'ro',ms = 6)    
xA1 = x[1] + fA; ax.plot(xA1,0,'ro',ms = 6)     
xB0 = x[3] - fB; ax.plot(xB0,0,'bo',ms = 6)    
xB1 = x[3] + fB; ax.plot(xB1,0,'bo',ms = 6)

# Output plane and focal plane
ax.plot(XB,YB0,'m',lw = 1) 
ax.plot(XB,YB1,'m',lw = 1) 
ax.plot(XB,YB2,'m',lw = 1) 
ax.plot(xA,yA,'ro',ms = 5) 
ax.plot(xB,yB,'bo',ms = 5) 

# Plane of lens
yLim = ax.get_ylim()
X = [x[1],x[1]]; Y = yLim; ax.plot(X,Y,'r',lw = 1) 
X = [x[3],x[3]]; Y = yLim; ax.plot(X,Y,'b',lw = 1)     

ax.set_xlim([X1,X2])     
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.set_title('Image point A: (%0.2f' %xA[0] + ', %0.2f) '  %yA[0] 
        + '   Image point B: ( %0.2f' %xB[0] + ', %0.2f) '  %yB[0] , fontsize = 10)
ax.grid()
fig1.tight_layout()

#%% OUTPUT TO CONSOLE 
print(' ')
print('3 ray paths through a two converging lens optical system')
print('Optical axis limits: X1 = %0.0f' % X1 + '   X2 = %0.0f' % X2)
print('focal lengths: fA = %0.2f' % fA + '   fB = %0.2f' % fB)
print('Translations: L01 = %0.1f ' % L01 + 
     '  L23 = %0.1f' % L23 + '  L45 = %0.1f ' % L45 )
print('Input (object): height y[0] = %0.2f' % y[0][0])
print('Input (object): elevations a0deg = %0.1f ' % a0deg[0] + 
     '  %0.1f' % a0deg[1] + '  %0.1f   deg ' %a0deg[2] )
print('Intermediate image lens A: xA = %0.3f' %xA[0] + '   height yA = %0.3f' %yA[0])
print('Output image lens B: xB = %0.3f' %xB[0] + '   height yB = %0.3f' %yB[0])
print('Transverse magnification:  mT = %0.3f' % mT[0])
print('Angular magnification:  mA = %0.3f' % mA)
q = xA[0] - x[2]; print('Lens A sA0 = %0.2f' % sA0 + '   sA1 = %0.2f' % q ) 
q0 = abs(xA[0] - x[4]); q1 = abs(xB[0] - x[4])
print('Lens B sB0 = %0.2f' % q0 + '   sB1 = %0.2f' % q1 )
tExe = time.time() - tStart
print('\nExecution time %0.0f s' %tExe)

#%% Save figures
# fig1.savefig('a1.png')   
     





