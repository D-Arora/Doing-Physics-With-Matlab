# -*- coding: utf-8 -*-
'''
ds2700.py      Nov 2025
DYNAMICAL SYSTEMS
  LORENZ EQUATIONS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS.ds2700.pdf
'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt, real, imag 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import find_peaks
from numpy.linalg import eig
import time

tStart = time.time()
plt.close('all')

#%% FUNCTIONS  Solve ODE for x,y and evaluate the Jacobian matrix   
def lorenz(t, state): 
    x, y, z = state
    dx = s*(y - x)
    dy = r*x - y -x*z
    dz = x*y - b*z
    return [dx, dy, dz]  

def JM(X,Y,Z,R):
    a11 = -s; a12 = s; a13 = 0
    a21 = R-Z; a22 = -1; a23 = -X
    a31 = Y; a32 = X; a33 = -b
    J = array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
    return J

#%%  INPUTS >>>
# System parameters
r = 28
s, b = 10, 8/3
N = 299
# Initial conditions  comment / uncomment
x0, y0, z0 = 0.2,0.2,0.2
col = [1,0,0]
col = [0,0,1]
# Phase space limits
L = 80

# time span
tS = 40
nT = 9999

#%% Solution ODE for x and y 
t = linspace(0,tS,nT)
u0 = [x0,y0,z0]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]
zS = sol[:,2]    

# Fixed points other than the Origin
xEp,xEm,yEp,yEm,zE = 0,0,0,0,0
if r > 1:
   ev1 = zeros(N); ev2 = ev1; ev3 = ev1 
   xEp = sqrt(b*(r-1)); xEm = - xEp
   yEp = xEp; yEm = xEm; zE = r-1
   J = JM(xEp,yEp,zE,r)
   ev, eV = eig(J)
   ev1 = ev[0]; ev2 = ev[1]; ev3 = ev[2]

# Console output
print('   ')
print('Bifurcation parameter r = %0.3f' %r)
print('fixed points: xE =  %0.3f' %xEp + '     yE = %0.3f' %yEp
      + '     zE = %0.3f' %zE)
print('              xE = %0.3f' %xEm + '     yE = %0.3f' %yEm
      + '    zE = %0.3f' %zE)
if r > 1:
   print('Eigenvalues:')
   print('lambda1 = %0.3f' %real(ev1) + ' + %0.3f j' %imag(ev1))
   print('lambda2 = %0.3f' %real(ev2) + ' + %0.3f j' %imag(ev2))
   print('lambda3 = %0.3f' %real(ev3) + ' + %0.3f j' %imag(ev3))
   
 
print('Initial conditions: x(0) = %0.3f' %x0 + '     y(0) = %0.3f' %y0
      + '     z(0) = %0.3f' %z0)   
print('Simulation time tS = %0.2f' %tS)

print('End points: xF = %0.3f' %xS[-1] + '     yF = %0.3f' %yS[-1]
      + '     zF = %0.3f' %zS[-1])



#%% FIG 1:  3D trajectory plot
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,4)
fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
ax.set_xlabel('x'); ax.set_ylabel('y');ax.set_zlabel('z')
ax.set_title('r = %0.3f' %r,fontsize = 10)

ax.plot3D(xS,yS,zS,lw = 1, color = col)
ax.plot3D(xS[0],yS[0],zS[0],'go',ms = 7)
ax.plot3D(xEp,yEp,zE,'ko',ms = 7)
ax.plot3D(xEm,yEm,zE,'ko',ms = 7)
ax.plot3D(0,0,0,'ko',ms = 7)

fig1.tight_layout()
fig1.savefig('a1.png')

#%% FIG 2: x,y,z trajectories
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,5)
fig2, ax = plt.subplots(nrows=3, ncols=1)

ax[0].set_xlabel('t'); ax[0].set_ylabel('x', rotation = 0)
ax[0].set_title('r = %0.2f' % r)
ax[0].grid()
ax[0].set_ylim(-max(xS),max(xS))
ax[0].plot(t,xS,lw = 2,color = col)
ax[0].plot([0,tS],[0,0],lw = 1,color = 'm')
ax[0].plot([0,tS],[xEp,xEp],lw = 1,color = 'm')
ax[0].plot([0,tS],[xEm,xEm],lw = 1,color = 'm')
 
ax[1].set_xlabel('t'); ax[1].set_ylabel('y', rotation = 0)
ax[1].grid()
ax[1].set_ylim(-max(yS),max(yS))
ax[1].plot(t,yS,lw = 2,color = col)
ax[1].plot([0,tS],[0,0],lw = 1,color = 'm')
ax[1].plot([0,tS],[yEm,yEm],lw = 1,color = 'm')
ax[1].plot([0,tS],[yEp,yEp],lw = 1,color = 'm')

ax[2].set_xlabel('t'); ax[2].set_ylabel('z', rotation = 0)
ax[2].grid()
ax[2].set_ylim(-max(zS),max(zS))
ax[2].plot(t,zS,lw = 2,color = col)
ax[2].plot([0,tS],[0,0],lw = 1,color = 'm')
ax[2].plot([0,tS],[zE,zE],lw = 1,color = 'm')

fig2.tight_layout()
fig2.savefig('a2.png')

#%% FIG 3:  x trajectory plot
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('t'); ax.set_ylabel('x', rotation = 0)
ax.set_title('r = %0.2f' % r)
ax.grid()
ax.set_ylim(-max(xS),max(xS))
ax.plot(t,xS,lw = 2,color = col)
ax.plot([0,tS],[0,0],lw = 1,color = 'm')
ax.plot([0,tS],[xEp,xEp],lw = 1,color = 'm')
ax.plot([0,tS],[xEm,xEm],lw = 1,color = 'm')
fig3.tight_layout()
fig3.savefig('a3.png')

#%% FIG 4: Phase portraits
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,6)
fig4, ax = plt.subplots(nrows=3, ncols=1)

ax[0].set_xlabel('x'); ax[0].set_ylabel('y', rotation = 0)
ax[0].set_title('r = %0.2f' % r)
ax[0].grid()
#ax[0].set_ylim(-max(xS),max(xS))
ax[0].plot(xS,yS,'b',lw = 0.5)
ax[0].plot(0,0,'ro',ms = 5)
ax[0].plot(xEp,yEp,'ro',ms = 5)
ax[0].plot(xEm,yEm,'ro', ms = 5)

ax[1].set_xlabel('x'); ax[1].set_ylabel('z', rotation = 0)
ax[1].grid()
#ax[0].set_ylim(-max(xS),max(xS))
ax[1].plot(xS,zS,'b',lw = 0.5)
ax[1].plot(0,0,'ro',ms = 5)
ax[1].plot(xEp,zE,'ro',ms = 5)
ax[1].plot(xEm,zE,'ro', ms = 5) 

ax[2].set_xlabel('y'); ax[2].set_ylabel('z', rotation = 0)
ax[2].grid()
#ax[0].set_ylim(-max(xS),max(xS))
ax[2].plot(yS,zS,'b',lw = 0.5)
ax[2].plot(0,0,'ro',ms = 5)
ax[2].plot(yEp,zE,'ro',ms = 5)
ax[2].plot(yEm,zE,'ro', ms = 5) 

fig4.tight_layout()
fig4.savefig('a4.png')

#%% fixed points C+ C-    r > 1
N = 999
R = linspace(1.1,40,N)
zE = R-1
xE1 = sqrt(b*(R-1)); yE1 = xE1 
xE2 = -xE1; yE2 = xE2 
zE = (R-1)*ones(N)

EV1 = zeros([N,3]); EV2 = zeros([N,3])
for c in range(N):
    J = JM(xE1[c],yE1[c],zE[c],R[c])
    ev, eV = eig(J)
    EV1[c,0]= real(ev[0]); EV1[c,1] = real(ev[1]); EV1[c,2] = real(ev[2])
    EV2[c,0]= imag(ev[0]); EV2[c,1] = imag(ev[1]); EV2[c,2] = imag(ev[2])

lambda1 = EV1[:,1]
rH = R[lambda1>0][0]

# FIG 5:
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r')
ax.set_ylabel(r'fixed points:  x$_E$   y$_E$   z$_E$',fontsize = 10)

ax.grid()
ax.plot(R,xE1,'r',lw = 2,label ='C$^+$')
ax.plot(R,zE,'b', lw = 2,label = '$z_E$')
ax.plot(R,xE2,'m',lw = 2,label = 'C$^-$')
ax.legend(ncols=3)
fig5.tight_layout()
fig5.savefig('a5.png')

# FIG 6:
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig6, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r')
ax.set_ylabel('$\lambda$',rotation = 0)
ax.set_title('Fixed points: C$^+$   C$^-$',fontsize = 12)
ax.text(26,-4,'r$_H$ = %0.2f' %rH)
ax.grid()
ax.plot(R,EV1[:,0],'b',lw = 2)
ax.plot(R,EV1[:,1],'r', lw = 2)
ax.plot(R,EV1[:,2],'m',lw = 2)
ax.plot([rH,rH],[-15,1],'k',lw=1)
fig6.tight_layout()
fig6.savefig('a6.png')


#%% Stability at the Origin
EV = zeros([N,3])
R = linspace(0.1,5,N)
for c in range(N):
    J = JM(0,0,0,R[c])
    ev, eV = eig(J)
    EV[c,0]= ev[0]; EV[c,1] = ev[1]; EV[c,2] = ev[2]

# FIG 7:    
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig7, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r')
ax.set_ylabel('$\lambda$', rotation = 0)
ax.set_title('Fixed point: Origin (0,0,0)',fontsize = 10)
ax.grid()
ax.plot(R,EV[:,0],'b',lw = 2)
ax.plot(R,EV[:,1],'r',lw = 2)
ax.plot(R,EV[:,2],'k',lw = 2)
fig7.tight_layout()
fig7.savefig('a7.png')


#%%
def lorenzR(t, state): 
    x, y, z = state
    dx = s*(y - x)
    dy = p*x - y - x*z
    dz = x*y - b*z
    return [dx, dy, dz]  

col = [1,0,0]
col = [0,0,1]
N = 999
tS = 20; nT = 999; t = linspace(0,tS,nT)
#R = linspace(0.1, 40, N)
R = linspace(0.1, 40, N)
u0 = [2.00,2,20]
xF = zeros(N); yF = zeros(N); zF = zeros(N)
xEp = zeros(N)

for c in range(N):
    p = R[c]
    sol = odeint(lorenzR, u0, t, tfirst=True)
    xS = sol[:,0]     
    yS = sol[:,1]
    zS = sol[:,2]    
    xF[c] = xS[-1]; yF[c] = zF[c] = zS[-1]
    if R[c]>1: xEp[c]= sqrt(b*(R[c]-1))

# FIG 8
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig8, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r')
ax.set_ylabel('xF')
ax.grid()
ax.plot(R,xEp,'m',lw = 0.5)
ax.plot(R,-xEp,'m', lw =0.5)
#ax.plot([0,R[-1]],[0,0],'m', lw =0.5)
ax.plot(R,xF,'o',color = col, ms = 2)

fig8.tight_layout()
fig8.savefig('a8.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


