# -*- coding: utf-8 -*-
"""
Nov 2026         op200.py

OPTICS
    FERMAT'S PRINCIPLE
    LAW OF REFLECTION and LAW OF REFRACTION Snell's Law'
    

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/op200.pdf

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time
from scipy.integrate import solve_ivp, odeint, simpson
from numpy import pi, linspace, zeros, array, ones, sqrt, cos, sin, exp, real, imag

tStart = time.time()
plt.close('all')

#%% FUNCTION
def dist(x1,y1,x2,y2):
    L = sqrt((x1-x2)**2 + (y1-y2)**2)
    return L

#%% REFLECTION: Find yB
N = 9999
n1 = 1
xA = 0; yA = -1; xB = 1; xC = 0; yC = 0.8

y = linspace(yA,yC,N)
opl_AB = n1*dist(xA,yA,xB,y)
opl_BC = n1*dist(xC,yC,xB,y)
opl_ABC = opl_AB + opl_BC
index = np.where(opl_ABC == min(opl_ABC))
yB = y[index]; yB = yB[0]
opl_ABC_min = opl_ABC[index]
theta1 = np.degrees(np.arctan((yB-yA)/xB))
theta2 = np.degrees(np.arctan((yC-yB)/xB))

#%% CONSOLE DISPLAY
print('REFLECTION')
print('  Angle of incidence,  theta1 = %0.2f deg' %theta1 )
print('  Angle of reflection, theta2 = %0.2f deg' %theta2 )

#%% REFLECTION  OPL
txtX = 'y' ; txtY = 'opl'
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
xP = y; yP = opl_ABC; ax.plot(xP,yP,'b',lw = 2)
xP = yB; yP = opl_ABC_min; ax.plot(xP,yP,'or',ms = 5)
ax.set_title('y$_B$ = %0.2f' %yB)

fig1.tight_layout()
fig1.savefig('a1.png')

#%%  REFLECTION RAY DIAGRAM
txtX = 'x' ; txtY = 'y'
plt.rcParams['font.size'] = 12
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (4,4))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()

xP = [xA,xB]; yP = [yA,yB]; ax.plot(xP,yP,'b',lw = 2)
xP = [xB,xC]; yP = [yB,yC]; ax.plot(xP,yP,'r',lw = 2)
xP = [xA,xB]; yP = [yB,yB]; ax.plot(xP,yP,'k',lw = 1)

xP = [xA,xB]; yP = [yA,y[1200]]; ax.plot(xP,yP,'m',lw = 1)
xP = [xB,xC]; yP = [y[1200],yC]; ax.plot(xP,yP,'m',lw = 1)

xP = [xA,xB]; yP = [yA,y[8000]]; ax.plot(xP,yP,'m',lw = 1)
xP = [xB,xC]; yP = [y[8000],yC]; ax.plot(xP,yP,'m',lw = 1)
xP = [xB,xB]; yP = [-1,1]; ax.plot(xP,yP,'k',lw = 5)

ax.text(-0.32,yB-0.2,r'$\theta_1$ = %0.1f$^o$' %theta1,color = 'b' )
ax.text(-0.32,yB+0.2,r'$\theta_2$ = %0.1f$^o$' %theta2,color = 'r' )
ax.text(-0.15,yA,'A')
ax.text(xB+0.05,yB,'B')
ax.text(-0.15,yC,'C')

fig2.tight_layout()
ax.axis('equal')  
ax.set_yticks(np.arange(-1,1.1,0.5))
fig2.savefig('a2.png')


########################################################
# REFRACTION
########################################################
# Find yB
yD = 0.8    # <<< input value  0 < yD < 0.8

xD = 2
n1 = 1; n3 = 1.5 
# CALCULATIONS
y = linspace(yA,yD,N)
opl_AB = n1*dist(xA,yA,xB,y)
opl_BD = n3*dist(xD,yD,xB,y)
opl_ABD = opl_AB + opl_BD
index = np.where(opl_ABD == min(opl_ABD))
yB = y[index]; yB = yB[0]
opl_ABD_min = opl_ABD[index]

theta1 = np.degrees(np.arctan((yB-yA)/xB))
theta3 = np.degrees(np.arctan((yD-yB)/(xD-xB)))

LHS = n1*sin(np.radians(theta1))
RHS = n3*sin(np.radians(theta3))

print('  ')
print('RERACTION')
print('  Angle of incidence,  theta1 = %0.2f deg' %theta1 )
print('  Angle of refraction, theta2 = %0.2f deg' %theta2 )
print('  LHS = %0.2f ' %LHS )
print('  LHS = %0.2f ' %RHS )


# REFRACTION  OPL
txtX = 'y' ; txtY = 'opl'
plt.rcParams['font.size'] = 12
fig3, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.5))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.set_title('y$_B$ = %0.3f' %yB)
ax.grid()
xP = y; yP = opl_ABD; ax.plot(xP,yP,'b',lw = 2)
xP = yB; yP = opl_ABD_min; ax.plot(xP,yP,'or',ms = 5)
fig3.tight_layout()
fig3.savefig('a3.png')

#%%  REFRACTION RAY DIAGRAM
txtX = 'x' ; txtY = 'y'
plt.rcParams['font.size'] = 12
fig4, ax = plt.subplots(nrows=1, ncols=1,figsize = (4,4))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()

xP = [xA,xB]; yP = [yA,yB]; ax.plot(xP,yP,'b',lw = 2)
xP = [xB,xD]; yP = [yB,yD]; ax.plot(xP,yP,'r',lw = 2)
xP = [xA,xD]; yP = [yB,yB]; ax.plot(xP,yP,'k',lw = 1)

xP = [xA,xB]; yP = [yA,y[1200]]; ax.plot(xP,yP,'m',lw = 1)
xP = [xB,xD]; yP = [y[1200],yD]; ax.plot(xP,yP,'m',lw = 1)
xP = [xA,xB]; yP = [yA,y[8000]]; ax.plot(xP,yP,'m',lw = 1)
xP = [xB,xD]; yP = [y[8000],yD]; ax.plot(xP,yP,'m',lw = 1)

ax.text(0.12,yB-0.2,r'$\theta_1$ = %0.1f$^o$' %theta1,color = 'b' )
ax.text(1.4,yB+0.06,r'$\theta_3$ = %0.1f$^o$' %theta2,color = 'r' )

ax.text(0.1,yA,'A')
ax.text(xB,yB-0.15,'B')
ax.text(xD-0.1,yD-0.2,'D')

ax.axis('equal')  
ax.set_yticks(np.arange(-1,1.1,0.5))
ax.set_xlim([0,2]); ax.set_ylim([-1,1])
ax.axvspan(1,2, color="orange", alpha=0.2)  # Vertical band from x=1 to x=2

fig4.tight_layout()
fig4.savefig('a4.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)


