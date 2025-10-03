# -*- coding: utf-8 -*-

'''
ds1510.py      OCt 2025
NONLINEAR [2D] DYNAMICAL SYSTEMS
SUBCRITICAL HOPF BIFURCATION

Ian Cooper: matlabvisualphysics@gmail.com

   https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
   https://d-arora.github.io/Doing-Physics-With-Matlab/pyDS1510.pdf

'''
# Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, sqrt 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')
tStart = time.time()

#%% FUNCTIONS  Solve ODE for x,y    x = R   y = theta 
def lorenz(t, state):    
    x, y = state
    dx = r*x + x**3 - x**5
    dy = 1
    return [dx, dy]  

#%%
# Bifurcation parameter  >>>
r = 0.25

# Input initial conditions  radius R and polar angle T (theta) >>>
R0,T0 = 0.1, -0.25*pi

col = [0,0,1]

#%% Solution ODE for x and y 
t1 = 0; t2 = 12; nT = 999
t = linspace(t1,t2,nT)
u0 = [R0,T0]
sol = odeint(lorenz, u0, t, tfirst=True)
R = sol[:,0]     
T = sol[:,1]            
xS = R*cos(T)
yS = R*sin(T)

# Steady-state radii
RP = ((1+(1+4*r)**0.5)/2)**0.5
RM = ((1-(1+4*r)**0.5)/2)**0.5

print(RP)
print(RM)


#%%  FIG 1: Phase portrait and trajectories
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,4)
fig1, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.3f' % r, fontsize = 14)
axes.set_xlim([-1.2, 1.2])
axes.set_ylim([-1.2, 1.2])
axes.grid()

axes.plot(xS[0], yS[0],'go',ms = 8)
axes.plot(xS, yS,'b',lw = 2)

A = linspace(0,2*pi,999)
xP = RM*cos(A); yP = RM*sin(A)
axes.plot(xP,yP,'k',lw = 3)
xP = RP*cos(A); yP = RP*sin(A)
axes.plot(xP,yP,'k',lw = 3)

axes.set_aspect('equal', 'box')
fig1.tight_layout()
fig1.savefig('a1.png')


#%% FIG 2:Phase Portrait streamplot  
x1 = -2; x2 = 2; nX = 21
xP = linspace(x1,x2,nX)
y1 = -2; y2 = 2; nY = 21
yP = linspace(y1,y2,nX)
xx,yy = np.meshgrid(xP,yP)
xxDot = r*xx + xx**3 - xx**5
yyDot = -yy

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,4)
fig2, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('x',color= 'black',fontsize = 12)
axes.set_ylabel('y',color = 'black',fontsize = 12)
axes.set_title('r = %2.3f' % r, fontsize = 14)
axes.set_xlim([-2, 2])
axes.set_ylim([-2, 2])

axes.plot(0,0,'ko',ms = 6)
axes.plot(xS[0], yS[0],'go',ms = 8)
axes.plot(xS, yS,color = col,lw = 2)

axes.streamplot(xx,yy,xxDot,yyDot, density = 1.2,
                linewidth = 1,color = 'k')

axes.set_aspect('equal', 'box')
fig2.tight_layout()
fig2.savefig('a2.png')


#%%  FIG 3: TIME EVOLUTION PLOTS  T r x y

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3.6)
fig3, axes = plt.subplots(nrows=2, ncols=2)
fig3.suptitle('r = %2.3f' % r, fontsize = 14)  

Rg = 0;Cg = 0   
axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Rg,Cg].set_ylabel(r'$\theta$',color = 'black',fontsize = 12)
# #axes[C].set_xlim([-5, 5])
# #axes[C].set_xticks(np.arange(-5,5.1,1))
#axes[Rg,Cg].set_ylim([0, 26])
# #axes[R,C].set_yticks(np.arange(-20,81,20))
axes[Rg,Cg].xaxis.grid()
axes[Rg,Cg].yaxis.grid()
axes[Rg,Cg].plot(t, T, color = col,lw = 2)

Rg = 0;Cg = 1   
axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Rg,Cg].set_ylabel('R',color = 'black',fontsize = 12)

axes[Rg,Cg].xaxis.grid()
axes[Rg,Cg].yaxis.grid()
axes[Rg,Cg].plot(t, R,color = col,lw = 2)

Rg = 1;Cg = 0   
axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Rg,Cg].set_ylabel('x',color = 'black',fontsize = 12)
axes[Rg,Cg].xaxis.grid()
axes[Rg,Cg].yaxis.grid()
axes[Rg,Cg].plot(t, xS,color = col,lw = 2)

Rg = 1;Cg = 1   
axes[Rg,Cg].set_xlabel('t',color= 'black',fontsize = 12)
axes[Rg,Cg].set_ylabel('y',color = 'black',fontsize = 12)

axes[Rg,Cg].xaxis.grid()
axes[Rg,Cg].yaxis.grid()
axes[Rg,Cg].plot(t, yS,color = col,lw = 2)

fig3.tight_layout()
fig3.savefig('a3.png')

#%%   FIG 4:  Rdot as a function of r >>>
R1 = 0; R2 = 1; N = 599
rp = -0.15

Rp = linspace(R1,R2,N)
Rdot = rp*Rp + Rp**3 - Rp**5

RP = ((1+(1+4*rp)**0.5)/2)**0.5
RM = ((1-(1+4*rp)**0.5)/2)**0.5

fRP = r + 3*RP**2 -5*RP**4
fRM = r + 3*RM**2 -5*RM**4

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (4,3)
fig4, axes = plt.subplots(nrows=1, ncols=1)

axes.set_xlabel('R',color= 'black',fontsize = 12)
axes.set_ylabel('R$_{dot}$',color = 'black',fontsize = 12)
axes.set_title('r = %2.3f' % rp, fontsize = 14)
axes.grid()

axes.plot(0,0,'bo',ms = 6)
axes.plot(RM,0,'ro',ms = 6)
axes.plot(RP,0,'bo',ms = 6)
axes.plot(Rp,Rdot,'k',lw = 2)

fig4.tight_layout()
fig4.savefig('a4.png')


#%% Figure 5:  R vs  Rdot for different r values
def RDOT(rp):
    Rdot = rp*Rp + Rp**3 - Rp**5
    return Rdot

R1 = 0; R2 = 1.2; N = 599; Rp = linspace(R1,R2,N)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('R'); ax.set_ylabel('R$_{dot}$')
ax.grid()

rP = -0.5; Rdot = RDOT(rP); txt = str(rP)
ax.plot(Rp,Rdot,lw = 2, label = txt)

rP = -0.15; Rdot = RDOT(rP); txt = str(rP)
Re1 = sqrt((1+sqrt(1+4*rP))/2)
Re2 = sqrt((1-sqrt(1+4*rP))/2)
ax.plot(Rp,Rdot,lw = 2, label = txt)
ax.plot(Re1,0,'ko',ms = 3)
ax.plot(Re2,0,'ko',ms = 3)
ax.plot(0,0,'ko',ms = 3)

rP = 0; Rdot = RDOT(rP); txt = str(rP)
ax.plot(Rp,Rdot,lw = 2, label = txt)
rP = 0.25; Rdot = RDOT(rP); txt = str(rP)
ax.plot(Rp,Rdot,lw = 2, label = txt)
rP = 0.5; Rdot = RDOT(rP); txt = str(rP)
ax.plot(Rp,Rdot,lw = 2, label = txt)
ax.plot([0,1.2],[0,0],'k',lw = 1)

ax.legend(fontsize = 10)
fig5.tight_layout()
fig5.savefig('a5.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)
