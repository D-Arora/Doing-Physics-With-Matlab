# -*- coding: utf-8 -*-
"""
ds25L9.py
Sep 25

DYNAMICAL SYSTEMS:
    IMPERFECT Pitchfork Bifurcations

# Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/dsDocs/ds25L9.pdf


"""

# LIBRARIES  ================================================================
import numpy as np
from numpy import pi, sin, cos, linspace, sqrt, array, zeros, real, imag
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time

plt.close('all')
tStart = time.time()

#%% CELL 1 Fig 1:  fixed points  r = 27 > 0
# INPUT r value  r = 25 / x range 
N = 999
# INPUT r values >>>
r = 27

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('r = %0.0f' %r)
ax.grid()

if r > 0:
   x1 = -6.5; x2 = 6.5
   h = array([-25,25,-75,75])
   x = linspace(x1,x2,N)
   y2 = r*x - x**3 
   ax.set_xticks([-6,-3,0,3,6])
   ax.plot(x,y2,'k',lw = 2)
  
   xP = [x1,x2]
   y1 = [h[0],h[0]]; ax.plot(xP,y1,'b',lw = 2)
   y1 = [h[1],h[1]]; ax.plot(xP,y1,'b',lw = 2)
   y1 = [h[2],h[2]]; ax.plot(xP,y1,'m',lw = 2)
   y1 = [h[3],h[3]]; ax.plot(xP,y1,'m',lw = 2)

# Turning points: Peak and trough in cubic equation 
   xC = sqrt(r/3); Y1 = r*xC - xC**3
   yP = [Y1,Y1]; ax.plot(xP,yP,'r',lw = 2) 
   yP = [-Y1,-Y1]; ax.plot(xP,yP,'r',lw = 2) 
   hC = 2*(r/3)**1.5
   xC = sqrt(r/3)
   print('')
   print('CELL 1  ')
   print('rC = %0.2f' %r + '   hC = %0.5f ' %-hC + '   xC = %0.2f' %xC)
   print('rC = %0.2f' %r + '   hC = %0.5f ' %hC + '   xC = %0.2f' %-xC)

if r <= 0:
   x1 = -4; x2 = 4
   h = array([-150,-50,50,150])
   x = linspace(x1,x2,N)
   y2 = r*x - x**3 
   #ax.set_xticks([-6,-3,0,3,6])
   ax.plot(x,y2,'k',lw = 2) 
   xP = [x1,x2]
   y1 = [h[0],h[0]]; ax.plot(xP,y1,'b',lw = 2)
   y1 = [h[1],h[1]]; ax.plot(xP,y1,'b',lw = 2)
   y1 = [h[2],h[2]]; ax.plot(xP,y1,'b',lw = 2)
   y1 = [h[3],h[3]]; ax.plot(xP,y1,'b',lw = 2)
   
fig1.tight_layout()
fig1.savefig('a1.png')


#%% CELL 2   Fig 2:  fixed points  r = -25 <= 0
# INPUT r value  r = 25 / x range / h values
r = 8.77205; x1 = -4; x2 = 4
h = array([-250,-100,100,250])
xP = linspace(x1,x2,N)
y2 = r*x - x**3

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('r = %0.0f' %r)
ax.grid()

ax.plot(xP,y2,'k',lw = 2)

xP = [x1,x2]
y1 = [h[0],h[0]]; ax.plot(xP,y1,'b',lw = 2)
y1 = [h[1],h[1]]; ax.plot(xP,y1,'b',lw = 2)
y1 = [h[2],h[2]]; ax.plot(xP,y1,'b',lw = 2)
y1 = [h[3],h[3]]; ax.plot(xP,y1,'b',lw = 2)

fig2.tight_layout()
fig2.savefig('a2.png')


#%% CELL 3   Fig 3:    (h,r) diagram: number of fixed points - critical values for h
x1 = 0; x2 = 4
r = linspace(x1,x2,N)
hC = real(2*(r/3)**1.5)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r')
ax.set_ylabel('h$_C$')
#ax.grid()
ax.set_xlim([-2,4])
ax.plot(r,hC,'b',lw = 2)
ax.plot(r,-hC,'b',lw = 2)
ax.fill_between(r, hC, -hC, color='lightblue', alpha=0.5, label='Filled Area')
ax.text(-1,0,'1',color = 'k',fontsize = 14)
ax.text(3,0,'3',color = 'k',fontsize = 14)
ax.text(2,1,'2',color = 'b',fontsize = 14,bbox=dict(facecolor='yellow', alpha=1))
fig3.tight_layout()
fig3.savefig('a3.png')

#%%  CELL 4   Fig 4:  Bifurcation diagrams    xe vs r
# Zeros of xDot 
h = -10
N = 399
r = linspace(-25,50,N)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r')
ax.set_ylabel('x$_e$')
ax.set_title(' h = %0.0f' % h)
ax.grid()

for c in range(N):
    coeff = [-1,0,r[c],h]     
    Z = np.roots(coeff)
    e0 = Z[0]; e1 = Z[1]; e2 = Z[2]
    f0 = r[c] - 3*e0**2; f1 = r[c] - 3*e1**2; f2 = r[c] - 3*e2**2
    col = [0,0,1]
    if imag(e0) == 0:
       if f0 > 0: col = [1,0,0]
       ax.plot(r[c],e0,'o',color = col,ms = 1)
    col = [0,0,1]
    if imag(e1) == 0:
       if f1 > 0: col = [1,0,0]
       ax.plot(r[c],e1, 'o',color = col,ms = 1)
    col = [0,0,1]
    if imag(e2) == 0:
       if f2 > 0: col = [1,0,0]
       ax.plot(r[c],e2,'o', color = col,ms = 1)
    
fig4.tight_layout()
fig4.savefig('a4.png')

#%%  CELL 5   Fig 5:  Bifurcation diagrams    xe vs h  
# Zeros of xDot 
r = -25 
N = 399
h = linspace(-200,200,N)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('h')
ax.set_ylabel('x$_e$')
ax.set_title(' r = %0.0f' % r)
ax.grid()

for c in range(N):
    coeff = [-1,0,r,h[c]]     
    Z = np.roots(coeff)
    e0 = Z[0]; e1 = Z[1]; e2 = Z[2]
    f0 = r - 3*e0**2; f1 = r - 3*e1**2; f2 = r - 3*e2**2
    col = [0,0,1]
    if imag(e0) == 0:
       if f0 > 0: col = [1,0,0]
       ax.plot(h[c],e0,'o',color = col,ms = 1)
     
    col = [0,0,1]
    if imag(e1) == 0:
       if f1 > 0: col = [1,0,0]
       ax.plot(h[c],e1, 'o',color = col,ms = 1)
     
    col = [0,0,1]
    if imag(e2) == 0:
      if f2 > 0: col = [1,0,0]
      ax.plot(h[c],e2,'o', color = col,ms = 1)
     
fig5.tight_layout()
fig5.savefig('a5.png')

#%% CELL 6   Fig. 6  time evolxution of the flow along a line
#  SOLVE ODE for x
def lorenz(t, state):    
    xS = state
    dx = h + r*xS - xS**3
    return dx  

# Input h, r, x(0) = 0, t2
h = -10
r = 8.77
x0 = 5
t2 = 50

# Fixed points
rC = 3*(abs(h/2))**(2/3)
coeff = [-1,0,r,h]     
e = np.roots(coeff)
c = 0; z0 = 99; z1 = 99; z2 = 99
if imag(e[0]) == 0: z0 = real(e[0]); c = c+1 
if imag(e[1]) == 0: z1 = real(e[1]); c = c+1 
if imag(e[2]) == 0: z2 = real(e[2]); c = c+1 

N = 9999; t1 = 0; t = linspace(t1,t2,N)

sol = odeint(lorenz, x0, t, tfirst=True)
xS = sol[:,0] 

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig6, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_title('h = %0.1f' %h + '    r = %0.3f ' %r)
ax.grid()

ax.plot(t, xS, 'k',lw = 2)

fig6.tight_layout()
fig6.savefig('a6.png')

# Console output
txt0 = '   stable'; txt1 = '   stable'; txt2 = '   stable'
S0 = r - 3*z0**2; 
S1 = r - 3*z1**2
S2 = r - 3*z2**2
if S0 > 0: txt0 ='   unstable'
if S1 > 0: txt1 ='   unstable'
if S2 > 0: txt2 ='   unstable'

print('  ')
print('CELL 6')
print('h = %0.1f' %h + '    r = %0.5f ' %r)
print('cusp point: hC = %0.2f' %h + '   rC =% 0.5f' %rC)
print('x0 = %0.2f' %x0)
print('xEND = %0.3f' %xS[-1])
print('Number of fixed points = %0.0f' %c)
if z0 < 99: print('xe = %0.3f' %z0 + txt0)
if z1 < 99: print('xe = %0.3f' %z1 + txt1)
if z2 < 99: print('xe = %0.3f' %z2 + txt2)


#%% CELL 7   Fig. 7   xDot vs x
# INPUT h and r >>>
h = 25
r = 27
N = 999
x1 = -8; x2 = 8; x = linspace(x1,x2,N) 
xDot = h + r*x - x**3
hC = 2*(r/3)**1.5
# Fixed points
coeff = [-1,0,r,h]     
e = np.roots(coeff)
c = 0; z0 = 99; z1 = 99; z2 = 99
if imag(e[0]) == 0: z0 = real(e[0]); c = c+1 
if imag(e[1]) == 0: z1 = real(e[1]); c = c+1 
if imag(e[2]) == 0: z2 = real(e[2]); c = c+1 

txt0 = '   stable'; txt1 = '   stable'; txt2 = '   stable'
S0 = r - 3*z0**2; 
S1 = r - 3*z1**2
S2 = r - 3*z2**2
if S0 > 0: txt0 ='   unstable'
if S1 > 0: txt1 ='   unstable'
if S2 > 0: txt2 ='   unstable'

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig7, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x')
ax.set_ylabel('x$_{dot}$')
ax.set_title('h = %0.1f' %h + '    r = %0.1f ' %r +
             '   h$_C$ = %0.2f' %hC)
ax.set_xticks(np.arange(x1,x2+1,2))
ax.grid()

ax.plot(x, xDot, 'k',lw = 2)
ax.plot([x1,x2], [0,0], 'r',lw = 2)

# Console output
print('  ')
print('CELL 7')
print('h = %0.1f' %h + '    r = %0.5f ' %r)
print('Number of fixed points = %0.0f' %c)
if z0 < 99: print('xe = %0.3f' %z0 + txt0)
if z1 < 99: print('xe = %0.3f' %z1 + txt1)
if z2 < 99: print('xe = %0.3f' %z2 + txt2)

fig7.tight_layout()
fig7.savefig('a7.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)



