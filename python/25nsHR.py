# -*- coding: utf-8 -*-
'''
AA\ns\pyCode\ns25_HR.py      July 2025
# COMPLEX SYSTEMS: NEUROSCIENCE
# HINDMARSH-ROSE MODEL BURSTING NEURONS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/ns25_HR.pdf

https://pmc.ncbi.nlm.nih.gov/articles/PMC9871131/

https://scicomp.stackexchange.com/questions/36013/numerical-computation-of-lyapunov-exponent


'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, real, imag, sqrt 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import random
from scipy.optimize import fsolve
from sympy import symbols, Eq, solve

tStart = time.time()

plt.close('all')


#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x, y, z = state
    I = 0
    if t > T1: I = I0
    if t > T2: I = 0
    dx = y + 3*x**2 - x**3 - z + I
    dy = 1 - 5*x**2 - y
    dz = r*(s*(x-xc) - z)
    return [dx, dy, dz]  


#%% INPUTS >>>
# Initial conditions
x1,y1,z1 = 0.618, -1, 0

# Adpation current
r,s = 0.001, 4  
# Time span for solving ODE
t1 = 0; t2 = 500; nT = 9999     
t = linspace(t1,t2,nT)         # time span for solving ODEs
# External current stimulus: off - on times
I0 = 0
T1 = 0; T2 = t2
Iext = zeros(nT); Iext[t>T1] = I0; Iext[t>T2] = 0     
        

#%%   CRITICAL POINTS  xC and yC  z = 0
# Intersection of x and y nullclines: cubic equation --> 3 roots  
c1 = 1; c2 = 2; c3 = 0; c4 = -1 - I0
coeff = [c1,c2,c3,c4]
xC = np.roots(coeff)
yC = 1 - 5*xC**2  
zC = zeros(3)  
# x, y, z = symbols('x y z')
# SS = solve([y + 3*x**2 - x**3,1 - 5*x**2 - y], x,y)
# SS0 = SS[0]
# SS1 = SS[1]
# SS2 = SS[2]
# xC0 = float(SS0[0]); yC0 = float(SS0[1])
# xC1 = float(SS1[0]); yC1 = float(SS1[1])
# xC2 = float(SS2[0]); yC2 = float(SS2[1])
# xC = array([xC0,xC1,xC2])
# yC = array([yC0,yC1,yC2])


#%% Solve ODEs 
xc = -0.5*(1+sqrt(5)) 
u0 = [xc,y1,z1]
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1] 
zS = sol[:,2]  


#%%  PHASE PORTRAIT  x vs y   /   nullclines      I0 = constant  z = 0
X = linspace(-2,2,16); Y = linspace(-20,4,16)
xx,yy = np.meshgrid(X,Y)
xxDot = yy + 3*xx**2 - xx**3 + I0 
yyDot = 1 - 5**xx**2 - yy

XXDot = xxDot/(sqrt(xxDot**2 + yyDot**2))
YYDot = yyDot/(sqrt(xxDot**2 + yyDot**2))

# Nullclines   I0 = constant    z = 0
xN = linspace(-2,2,999)
yNx = -3*xN**2 + xN**3 - I0
yNy = 1 - 5*xN**2


#%% Console output
print('Model parameters: r = %5.4f ' % r, 's = %5.4f' % s, 'I0 = %5.4f ' % I0)
print('Initial values: x0 = %5.4f ' % x1, 'y0 = %5.4f' % y1, 'z0 = %5.4f' % z1)
print('Time span t2 = %0.0f ' % t2)
print('  ')
print('Critical values and Eigenvalues')

for c in range(3):
    flag1 = 0
    print('xC = %2.3f ' % real(xC[c]), '%2.3f' % imag(xC[c]),'j', \
       '   yC = %2.3f ' % real(yC[c]), '%2.3f' % imag(yC[c]),'j', \
       '   zC = %2.3f ' % real(zC[c]), '%2.3f' % imag(zC[c]),'j')
print('  ') 
       
#%% STABILIY of critical points
J12 = 1; J22 = -1
for c in range(3):
    J11 = -3*xC[c]**2 + 6*xC[c]
    J21 = -10*xC[c]
    J = np.array([ [J11,J12],[J21,J22] ]) 
    eV,eF = eig(J)      # eigenvalues and eigenfunctions
   
    flag0 = imag(eV[0]) == 0 and (real(eV[0]) < 0 and real(eV[1]) < 0)
    flag1 = imag(eV[0]) == 0 and (real(eV[0]) > 0 or real(eV[1]) > 0)
    flag2 = imag(eV[0]) == 0 and (real(eV[0]) == 0 or  real(eV[1]) == 0)
    flag3 = imag(eV[0]) != 0 and (real(eV[0]) > 0 or real(eV[1]) > 0)
    flag4 = imag(eV[0]) != 0 and (real(eV[0]) < 0 and real(eV[1]) < 0)
      
    if flag0 == 1:
           
           print('xC = %0.3f' %xC[c] + '  yC = %0.3f' %yC[c] + '  STABLE') 
           print('eV = ', np.round(eV,3))
    if flag1 == 1:
           print('xC = %0.3f' %xC[c] + '  yC = %0.3f' %yC[c] + '  UNSTABLE') 
           print('eV = ', np.round(eV,3))
    if flag2 == 1:
           print('xC = %0.3f' %xC[c] + '  yC = %0.3f' %yC[c] + '  inconclusive') 
           print('eV = ', np.round(eV,3))
    if flag3 == 1:
           print('xC = %0.3f' %xC[c] + '  yC = %0.3f' %yC[c] + '  unstable: osc growth')
           print('eV = ', np.round(eV,3))
    if flag4 == 1:
           print('xC = %0.3f' %xC[c] + '  yC = %0.3f' %yC[c] + '  stable: osc decay')
           print('eV = ', np.round(eV,3))
    print('  ')


#%% GRAPHICS
# FIGURE 1: Membrane potential x vs time t
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('time  t',color= 'black',fontsize = 12)
ax.set_ylabel('membrane potential  x',color = 'blue',fontsize = 12)
ax.grid()
ax.plot(t,xS,'b',lw = 2)
ax.tick_params(axis='y', colors='b')
ax1 = ax.twinx()
ax1.set_ylabel('I$_{ext}$', color = 'r',fontsize = 12)
ax1.tick_params(axis='y', colors='r')
ax1.plot(t,Iext,'r',lw = 2)
fig1.tight_layout()

#%% FIGURE 2: recovery var. y & adaption var. z vs time
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('time  t',color= 'black',fontsize = 12)
ax.set_ylabel('recovery var. y',color = 'blue',fontsize = 12)
ax.grid()
ax.plot(t,yS,'b',lw = 2)
ax.tick_params(axis='y', colors='b')
ax1 = ax.twinx()
ax1.set_ylabel('adaption var. z', color = 'r',fontsize = 12)
ax1.tick_params(axis='y', colors='r')
ax1.plot(t,zS,'r',lw = 2)
fig2.tight_layout()


#%% Figure 3: Phase Portrait: streamplot  
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('membrane pot. x',color= 'black',fontsize = 12)
ax.set_ylabel('recovery var. y',color = 'black',fontsize = 12)
ax.set_xlim([-2,2])
ax.set_ylim([-20,3])

ax.streamplot(xx,yy,xxDot,yyDot)
ax.plot(xN,yNx,'b',lw = 3)           #  nullcine xDot
ax.plot(xN,yNy,'r',lw = 3)           #  nullcine yDot

# Initial conditions and trajectories
u0 = [x1,y1,z1]
t1 = 0; t2 = 520
t = linspace(t1,t2,nT)
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1] 
zS = sol[:,2]      

ax.plot(xS,yS,'m',lw = 2.5)

# Critical points
for c in range(3):
    ax.plot(xC[c],yC[c],'ko',ms = 8)

ax.plot(xS[0],yS[0],'m+',ms = 10)

fig3.tight_layout()   

#%% Figure 4: Phase Portrait: quiver plot 
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
  
ax.set_xlabel('membrane pot. x',color= 'black',fontsize = 12)
ax.set_ylabel('recovery var. y',color = 'black',fontsize = 12)
ax.set_xlim([-2,2])
ax.set_ylim([-20,3])

ax.quiver(xx,yy,XXDot,YYDot)
ax.plot(xN,yNx,'b',lw = 3)           #  nullcine xDot
ax.plot(xN,yNy,'r',lw = 3)           #  nullcine yDot

# Initial conditions and trajectories
u0 = [x1,y1,z1]
t1 = 0; t2 = 520
t = linspace(t1,t2,nT)
sol = odeint(lorenz, u0, t, tfirst=True)
xS = sol[:,0]     
yS = sol[:,1] 
zS = sol[:,2]      

ax.plot(xS,yS,'m',lw = 2.5)

# Critical points
for c in range(3):
    ax.plot(xC[c],yC[c],'ko',ms = 8)

ax.plot(xS[0],yS[0],'m+',ms = 10)
 
fig4.tight_layout()   

#%% Figure 5: 3D tralectories  t x,y,z
fig5 = plt.figure(figsize = (5, 3))
ax = plt.axes(projection='3d')
ax.plot3D(xS,yS,zS,'b',lw = 2)
ax.plot3D(xS[-1],yS[-1],zS[-1],'ro')
ax.plot3D(xS[0],yS[0],zS[0],'bo')
ax.set_xlabel('x')
ax.set_ylabel(' y')
ax.set_zlabel('  z');
fig5.tight_layout()   

#%% Figure 6: phase plots
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,3)
fig6, ax = plt.subplots(nrows=1, ncols=3)
    
C = 0   
ax[C].set_xlabel('x',color = 'black',fontsize = 12)
ax[C].set_ylabel('y',color = 'black',fontsize = 12)
ax[C].xaxis.grid(); ax[C].yaxis.grid()
ax[C].plot(xS,yS,'b',lw = 2)
ax[C].plot(xS[0],yS[0],'go',ms = 6)

C = 1   
ax[C].set_xlabel('x',color = 'black',fontsize = 12)
ax[C].set_ylabel('z',color = 'black',fontsize = 12)
ax[C].xaxis.grid(); ax[C].yaxis.grid()
ax[C].plot(xS,zS,'b',lw = 2)
ax[C].plot(xS[0],zS[0],'go',ms = 6)

C = 2   
ax[C].set_xlabel('y',color = 'black',fontsize = 12)
ax[C].set_ylabel('z',color = 'black',fontsize = 12)
ax[C].xaxis.grid(); ax[C].yaxis.grid()
ax[C].plot(yS,zS,'b',lw = 2)
ax[C].plot(yS[0],zS[0],'go',ms = 6)

fig6.tight_layout()
    
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

#%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')



'''
#