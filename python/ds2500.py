# -*- coding: utf-8 -*-
'''
ds2500.py      Oct 2025
DYNAMICAL SYSTEMS
SAADLE NODE BIFURCATIONS IN PLANR SYSTEMS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS.ds2500.pdf



'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.linalg import eig
import time

tStart = time.time()
plt.close('all')

#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x, y = state
    dx = -a*x + y
    dy = x**2/(1+x**2) - b*y
    return [dx, dy]  

#%%  INPUTS >>>
a = 0.45
b = 1
L = 2
# Initial conditions xI yI (vI)
xI,yI = 0.16,0.80
# time span
tS = 150; nT = 99999
# Plot color
col = [1,0,0]

#%% Solution ODE for x and y 
t = linspace(0,tS,nT)
u0 = [xI,yI]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]    

#%% PHASE PORTRAIT streamplot
N = 11
x = linspace(0,L,N)
xx,yy = np.meshgrid(x,x)
xxDot = -a*xx + yy
yyDot = xx**2/(1+xx**2) - b*yy

#%% NULLCLINES
xNx = linspace(0,L,599)
yNx = a*xNx
xNy = linspace(0,L,599)
yNy = xNy**2/(b*(1+xNy**2))

#%% Fixed points, Jacobian, eigenvalues, stability
def JM(x):
    a11 = -a; a12 = 1; a21 = 2*x / (1+x**2)**2; a22 = -b
    J = array([[a11,a12],[a21,a22]])
    return J

# Critical value of bifurcation parameter ac
ac = 1/(2*b)
# Fixed points
xe = zeros(3); ye = zeros(3)
if a > ac:
   xe[0] = 0; ye[0] = 0
   J = JM(0)
   ev, eV = eig(J)
   print('Fixed point (xe,ye):  xe = %0.2f' %xe[0] + '   ye = %0.2f' %ye[0])
   print('Eigenvalues: ', ev)

if a == ac:
   xe[1] = 1/(2*a*b)
   ye[1] = a*xe[1]
   J = JM(xe[1])
   ev, eV = eig(J)
   print('Fixed point (xe,ye):  xe = %0.2f' %xe[1] + '   ye = %0.2f' %ye[1])
   print('Eigenvalues: %0.3f' %ev[0] + '   %0.3f' %ev[1])
   
if a < ac:
   xe[1] = (1 + sqrt(1-4*a**2*b**2)) / (2*a*b)
   ye[1] = a*xe[1]
   J = JM(xe[1])
   ev, eV = eig(J)
   print('Fixed point (xe,ye):  xe = %0.2f' %xe[1] + '   ye = %0.2f' %ye[1])
   print('Eigenvalues: %0.3f' %ev[0] + '   %0.3f' %ev[1])  
   
   xe[2] =(1 - sqrt(1-4*a**2*b**2)) / (2*a*b)
   ye[2] = a*xe[2]
   J = JM(xe[2])
   ev, eV = eig(J)
   print('Fixed point (xe,ye):  xe = %0.2f' %xe[2] + '   ye = %0.2f' %ye[2])
   print('Eigenvalues: %0.3f' %ev[0] + '   %0.3f' %ev[1])
   
#%%   FIG 1: t vs x   t vs y (1.00, 0.50) 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig1, ax = plt.subplots(nrows=2, ncols=1)
 
ax[0].set_xlabel('t'); ax[0].set_ylabel('x')
ax[0].set_title('a$_c$ = %0.2f' %ac + '   a = %0.2f' %a)
ax[0].grid()
ax[0].plot(t,xS,lw = 2,color = col)
 
ax[1].set_xlabel('t'); ax[1].set_ylabel('y')
ax[1].set_title('a$_c$ = %0.2f' %ac + '   a = %0.2f' %a)
ax[1].grid()
ax[1].plot(t,yS,lw = 2,color = col)

fig1.tight_layout()
fig1.savefig('a1.png')

#%% FIGURE 2: Phase Portrait  quiver
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig2, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x')
axes.set_ylabel('y', rotation = 0)
axes.grid()
axes.set_xlim([-0.03,2]);axes.set_ylim([-0.03,1]);

axes.quiver(xx,yy,xxDot,yyDot,linewidth = 1)

axes.plot(xe,ye,'mo',ms = 6)
axes.plot(xNx,yNx,'b')
axes.plot(xNy,yNy,'r')

axes.plot(xS,yS,'g',lw = 3)
axes.plot(xS[0],yS[0],'go',ms = 8)
axes.plot(xS[-1],yS[-1],'go',ms = 6)

fig2.tight_layout()
fig2.savefig('a2.png')

#%% FIGURE 3: Phase Portrait  streamplot
N = 99
x = linspace(0,2,N)
xx,yy = np.meshgrid(x,x)
xxDot = -a*xx + yy
yyDot = xx**2/(1+xx**2) - b*yy

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x')
axes.set_ylabel('y')
axes.grid()
axes.set_xlim([-0.03,2]);axes.set_ylim([-0.03,1])

axes.streamplot(xx,yy,xxDot,yyDot,linewidth = 1, density = 0.8, color='b')

axes.plot(xe,ye,'mo',ms = 6)
axes.plot(xNx,yNx,'b')
axes.plot(xNy,yNy,'r')

axes.plot(xS,yS,'g',lw = 3)
axes.plot(xS[0],yS[0],'go',ms = 8)
axes.plot(xS[-1],yS[-1],'go',ms = 6)

fig3.tight_layout()
fig3.savefig('a3.png')


#%% BIFURCATION DIAGRAM
aB = linspace(0.1,ac,599)
xP = (1 + sqrt(1-4*aB**2*b**2)) / (2*aB*b)
xM = (1 - sqrt(1-4*aB**2*b**2)) / (2*aB*b)
yP = aB*xP; yM = aB*xM

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=2)
C = 0
ax[C].set_xlabel('a')
ax[C].set_ylabel('x$_e$', rotation = 0)
ax[C].grid()
ax[C].set_ylim([-1,10]);
ax[C].plot(aB,xP,'b')
ax[C].plot(aB,xM,'r')
ax[C].plot([0,1],[0,0],'b')
C = 1
ax[C].set_xlabel('a')
ax[C].set_ylabel('y$_e$', rotation = 0)
ax[C].grid()
#ax[C].set_ylim([-1,10]);
ax[C].plot(aB,yP,'b',lw = 2)
ax[C].plot(aB,yM,'r')
ax[C].plot([0,1],[0,0],'b',lw = 2)

fig4.tight_layout()
fig4.savefig('a4.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


