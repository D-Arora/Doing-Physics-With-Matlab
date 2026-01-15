# -*- coding: utf-8 -*-
'''
ds2501.py      Oct 2025
DYNAMICAL SYSTEMS
PITCHFORK BIFURCATIONS IN PLANR SYSTEMS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS.ds2501.pdf
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
    dx = a*x + y + sin(x)
    dy = x - y
    return [dx, dy]  

#%%  INPUTS >>>
# bifurcation parameter
a = -1.6
# Initial conditions xI yI (vI)                
xI,yI = 1,1
# time span     
tS = 50; nT = 99999   
# dimensions of system           
L = 2 
# Plot color
col = [0,0,1]

ac = -2

#%% Solution ODE for x and y 
t = linspace(0,tS,nT)
u0 = [xI,yI]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]    

#%% NULLCLINES
xN = linspace(-L,L,599)
yNx = -(a*xN + sin(xN))
yNy = xN

#%% Fixed points, Jacobian, eigenvalues, stability
def JM(a,x):
    a11 = a + cos(x); a12 = 1; a21 = 1; a22 = -1
    J = array([[a11,a12],[a21,a22]])
    return J


#%%   FIG 1: t vs x   t vs y (1.00, 0.50) 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig1, ax = plt.subplots(nrows=2, ncols=1)
 
ax[0].set_xlabel('t'); ax[0].set_ylabel('x')
ax[0].set_title('a = %0.2f' % a + '    a$_C$ = %0.2f' %ac)
ax[0].grid()
ax[0].plot(t,xS,lw = 2,color = col)
 
ax[1].set_xlabel('t'); ax[1].set_ylabel('y')
ax[1].grid()
ax[1].plot(t,yS,lw = 2,color = col)

fig1.tight_layout()
fig1.savefig('a1.png')

#%% FIGURE 2: Phase Portrait  quiver
N = 11
x = linspace(-L,L,N)
xx,yy = np.meshgrid(x,x)
xxDot = a*xx + yy + sin(xx)
yyDot = xx - yy

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig2, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x')
axes.set_ylabel('y', rotation = 0)
axes.set_title('a = %0.2f' % a + '    a$_C$ = %0.2f' %ac)
axes.grid()
axes.set_xlim([-L,L]);axes.set_ylim([-L,L]);

axes.quiver(xx,yy,xxDot,yyDot,linewidth = 1)

axes.plot(0,0,'mo',ms = 6)
axes.plot(xN,yNx,'b')
axes.plot(xN,yNy,'r')

axes.plot(xS,yS,'g',lw = 3)
axes.plot(xS[0],yS[0],'go',ms = 8)
axes.plot(xS[-1],yS[-1],'go',ms = 6)

fig2.tight_layout()
fig2.savefig('a2.png')

#%% FIGURE 3: Phase Portrait  streamplot
N = 99
x = linspace(-L,L,N)
xx,yy = np.meshgrid(x,x)
xxDot = a*xx + yy + sin(xx)
yyDot = xx - yy

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x')
axes.set_ylabel('y',rotation = 0)
axes.set_title('a = %0.2f' % a + '    a$_C$ = %0.2f' %ac)
axes.grid()
axes.set_xlim([-L,L]);axes.set_ylim([-L,L])

axes.streamplot(xx,yy,xxDot,yyDot,linewidth = 1, density = 0.8, color='k')

axes.plot(0,0,'mo',ms = 6)
axes.plot(xN,yNx,'b')
axes.plot(xN,yNy,'r')

axes.plot(xS,yS,'g',lw = 3)
axes.plot(xS[0],yS[0],'go',ms = 8)
axes.plot(xS[-1],yS[-1],'go',ms = 6)

fig3.tight_layout()
fig3.savefig('a3.png')

#%% BIFURCATION DIAGRAM
def JB(a,x):
    a11 = a + cos(x); a12 = 1; a21 = 1; a22 = -1
    J = array([[a11,a12],[a21,a22]])
    return J
N = 99
a = linspace(-3,2,N)
ev0 = zeros(N); ev1 = zeros(N)
for c in range(N):
    J = JB(a[c],0)
    ev, eV = eig(J)
    ev0[c] = ev[0]; ev1[c] = ev[1]

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,2.5)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('a')
ax.set_ylabel(r'$\lambda$', rotation = 0)
ax.grid()

ax.plot(a,ev0,'b',lw = 2)
ax.plot(a,ev1,'r',lw = 2)

fig4.tight_layout()
fig4.savefig('a4.png')


#%% Find zeros 
num = 99999
q = -1.60
xE = linspace(-3,2,num)
fn = xE*(q+1) + sin(xE)
#Q = [0,0,0]
Q = np.zeros(3)
p = 0
for c in range(num-2):
    q = fn[c]*fn[c+1]
    if q <= 0:
       Q[p] = c
       p = int(p+1)     
Q = [int(x) for x in Q]     
print(xE[Q])


J = JM(-1.60,1.66)
ev, eV = eig(J)
print(np.round(ev,3))


#%%   Fixed points F = 0
q = -1.6
x = linspace(-2,2,999)
F = x*(q+1) + sin(x) 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,2.5)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('F(x)', rotation = 0)
ax.set_title('a = %0.2f' % q )
ax.grid()
ax.plot(x,F,'b',lw = 2)
fig5.tight_layout()
fig5.savefig('a5.png')




#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


