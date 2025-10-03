# -*- coding: utf-8 -*-

'''
ds1510.py      OCt 2025
NONLINEAR [2D] DYNAMICAL SYSTEMS
HOMOCLINIC BIFURCATIONS

Ian Cooper: matlabvisualphysics@gmail.com

   https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
   https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS1511.pdf

'''


# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')
tStart = time.time()

# FUNCTIONS  Solve ODE for x,y    x = R   y = theta 
def lorenz(t, state):    
    x, y = state
    dx = y
    dy = r*y + x - x**2 + x*y
    return [dx, dy]  

#%% INPUTS >>> r   x(0)   y(0)
r  = -0.90
x0 = 0.10
y0 = 0.1


#%% FIXED POINTS
xe0, ye0 = 0, 0
xe1, ye1 = -1, 0
xe2, ye2 = 1, 0

# Solution ODE for x and y 
t1 = 0; t2 = 30; nT = 999
t = linspace(t1,t2,nT)
u0 = [x0,y0]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1]            

#%%  Jacobian matrix and eigenvalues
rJ = r

J0 = np.array([[0,1],[1,rJ]])
ev0, J0 = eig(J0)
print(ev0)
JP = np.array([[0,1],[-1,rJ+1]])
evP, efP1 = eig(JP)
print(evP)
JM = np.array([[0,1],[3,rJ-1]])
evM, efPM = eig(JM)
print(evM)


#%% FIG 1:Phase Portrait streamplot  
x1 = -2; x2 = 2; nX = 299
xP = linspace(x1,x2,nX)
y1 = -2; y2 = 2; nY = 299
yP = linspace(y1,y2,nX)
xx,yy = np.meshgrid(xP,yP)
xxDot = yy
yyDot = r*yy + xx +xx*yy -xx**2

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig1, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.3f' % r, fontsize = 14)
axes.set_xlim([-2, 2])
axes.set_ylim([-2, 2])

axes.plot([x1,x2],[0,0],'m',lw=1)
axes.plot([xe1,xe1],[x1,x2],'m',lw=1)
axes.plot([xe2,xe2],[x1,x2],'m',lw=1)
axes.plot([0,0],[x1,x2],'m',lw=1)
# axes.plot(xS[0], yS[0],'go',ms = 8)
#axes.plot(xS, yS,'r',lw = 2)

axes.streamplot(xx,yy,xxDot,yyDot, density = 1.2,
                linewidth = 1,color = 'k')

axes.set_aspect('equal', 'box')
fig1.tight_layout()
fig1.savefig('a1.png')


#%% FIG 2: Phase portrait orbits

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3,3)
fig2, axes = plt.subplots(nrows=1, ncols=1)
  
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.3f   r$_C$ = 0865' % r, fontsize = 12)

axes.set_xlim([-2,2])

axes.set_ylim([-2,2])

axes.xaxis.grid()
axes.yaxis.grid()

axes.plot(xS[0], yS[0],'go',ms = 6)
axes.plot(xS, yS,'b',lw = 2)
axes.plot(0, 0,'ko',ms = 6)
axes.plot(1, 0,'ko',ms = 6)
axes.plot(-1, 0,'ko',ms = 6)

axes.set_aspect('equal', 'box')
fig2.tight_layout()

fig2.savefig('a2.png')

#%%  FIG 3: rajectories t vs x.y
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (3,3)
fig3, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('t',color= 'black',fontsize = 12)
axes.set_ylabel('x, y',color = 'black',fontsize = 12)
axes.set_title('r = %2.3f   r$_C$ = 0865' % r, fontsize = 12)
#axes.set_xlim([-1.2, 1.2])
#axes.set_ylim([-1.2, 1.2])
axes.grid()

axes.plot(t, xS,'b',lw = 2)
axes.plot(t, yS,'r',lw = 2)

# axes.set_aspect('equal', 'box')
fig3.tight_layout()
fig3.savefig('a3.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
