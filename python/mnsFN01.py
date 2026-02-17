# -*- coding: utf-8 -*-
'''

mnsFN01.py      FEB 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

FITZHUGH - NAGUMO MODEL

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mnsFN.pdf

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

#%%
def VDOT(t,v,u):
    Vdot  = v - v**3/3 - u + Iext
    return Vdot    
def UDOT(t,v,u):
    Udot  = (v + a - b*u)/c
    return Udot   


#%% INPUTS AND MODEL PARAMETERS
N = 9999
Iext = 1.4
col = [0,0,1]
# Fitzhugh-Nagum0 model parameters
a = 0.7; b = 0.8; c = 12.5
#a = 0.3; b = 1.4; c = 20
# Phase space plot dimensions: X axis [-3 3] / Y axis [-2 3]
Lx = [-3, 3]; Ly = [-2, 3]; nX = 16
# Vector field 
x1 = linspace(Lx[0],Lx[1],nX)
x2 = linspace(Ly[0],Ly[1],nX)
[xx, yy] = np.meshgrid(x1,x2)
f = (xx - xx**3/3 - yy + Iext)   
g = (1/c)*(xx + a - b*yy) 
fs = f/sqrt(f**2 + g**2)    # unit vectors
gs = g/sqrt(f**2 + g**2)

# Nullclines
vN = linspace(Lx[0],Lx[1],599)
uv = vN - vN**3/3 + Iext
uu = (vN + a)/b

# Fixed points (vE,uE)
c1= -b/3; c2 = 0; c3 = b-1; c4 = -a + b*Iext
VE = np.roots([c1,c2,c3,c4]) 

vE = np.real(VE[np.imag(VE) == 0])
uE = (vE + a)/b

for q in range(len(vE)):
    print('vE = %0.2f' %vE[q] + '   uE = %0.2f' % uE[q])


    
#%% Solve ODE with Runge-Kutta
N = 999
v = zeros(N); u = zeros(N)
# Initial condtions -2.8, -1.8
v[0] = -2.8
u[0] = -1.8

# Time span
tS = 200
t = linspace(0,tS,N)
h = t[2] - t[1]
       
#%% Runge-Kutta
for n in range(N-1):
    k1 = VDOT(t[n],v[n],u[n])
    g1 = UDOT(t[n],v[n],u[n])
    
    k2 = VDOT(t[n]+0.5*h,v[n]+k1*h/2,u[n]+g1*h/2)
    g2 = UDOT(t[n]+0.5*h,v[n]+k1*h/2,u[n]+g1*h/2)
    
    k3 = VDOT(t[n]+0.5*h,v[n]+k2*h/2,u[n]+g2*h/2)
    g3 = UDOT(t[n]+0.5*h,v[n]+k2*h/2,u[n]+g2*h/2)
    
    k4 = VDOT(t[n]+1.0*h,v[n]+k3*h,u[n]+g3*h)
    g4 = UDOT(t[n]+1.0*h,v[n]+k3*h,u[n]+g3*h)
    
    v[n+1] = v[n] + h*(k1 + 2*k2 + 2*k3 + k4)/6
    u[n+1] = u[n] + h*(g1 + 2*g2 + 2*g3 + g4)/6 


#%%   FIG 1: quiver plot  u vs v   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4.5,2.8)
fig1, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('v'); ax.set_ylabel('u')
ax.set_title('I$_{ext}$ = %0.3f' % Iext + '    v$_E$ = %0.2f' %vE[0]
          + '   u$_E$ = %0.2f' %uE[0], fontsize = 12)
ax.grid()
ax.set_ylim(Ly)

ax.quiver(xx,yy, fs,gs, color = [0,0,0])
ax.plot(vN,uv,'r',lw = 2)
ax.plot(vN,uu,'m',lw = 2)
ax.plot(v,u,'b',lw = 2)
ax.plot(v[0],u[0],'go',ms = 8)
ax.plot(vE,uE,'ko', ms = 8)

fig1.tight_layout()
fig1.savefig('a1.png')

#%%   FIG 2:  streamplot u vs v   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4.5,2.8)
fig2, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('v'); ax.set_ylabel('u')
ax.set_title('I$_{ext}$ = %0.3f' % Iext + '    v$_E$ = %0.2f' %vE[0]
          + '   u$_E$ = %0.2f' %uE[0], fontsize = 12)
ax.grid()
ax.set_ylim(Ly)
ax.streamplot(xx,yy, fs,gs, density=0.8, linewidth=1,color = [0,0,0])
ax.plot(vN,uv,'r',lw = 2)
ax.plot(vN,uu,'m',lw = 2)
ax.plot(v,u,'b',lw = 2)
ax.plot(v[0],u[0],'go',ms = 8)
ax.plot(vE,uE,'ko', ms = 8)

fig2.tight_layout()
fig2.savefig('a2.png')

#%%   FIG 3: t vs v  u   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.8)
fig3, ax = plt.subplots(nrows=1, ncols=1)
 
ax.set_xlabel('t  '); ax.set_ylabel('v & u')
ax.set_title('I$_{ext}$ = %0.3f' % Iext + '    v$_E$ = %0.2f' %vE[0]
          + '   u$_E$ = %0.2f' %uE[0], fontsize = 12)
ax.grid()
ax.plot(t,v,'b', lw = 2, label = 'v')
ax.plot(t,u,'r', lw = 2, label = 'u')
ax.legend()
fig3.tight_layout()
fig3.savefig('a3.png')

#%%   FIG 4: t vs IEXT   
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (5,2)
# fig4, ax = plt.subplots(nrows=1, ncols=1)
 
# ax.set_xlabel('t  [ ms ]'); ax.set_ylabel(r'I$_{ext}$  [ pA ]')
# ax.set_title('I$_{ext}$ = %0.0f  pA' %Iext)
# ax.grid()
# ax.set_ylim([0,1.1*max(Iext)])
# ax.plot(t,Iext,'r',lw = 2)

# fig4.tight_layout()
# fig4.savefig('a4.png')

# #%%   FIG 3: t vs u   
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (6,3)
# fig3, ax = plt.subplots(nrows=1, ncols=1)
 
# ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('u  [ mV ]')

# ax.set_title('C = %0.0f' %C+ '    k = %0.2f' %k +
#              '    v$_r$ = %0.0f' %vr + '    v$_t$ = %0.0f' %vt + 
#              '\n a = %0.2f' %a + '   b = %0.1f' %b +
#             '   c = %0.0f' %c + '    d = %0.0f' %d +
#             '   I$_0$ = %0.1f' %I0 , fontsize = 11)
# ax.grid()

# ax.plot(t,u,'b', lw = 2)

# fig3.tight_layout()
# fig3.savefig('a3.png')


# #%%   FIG 4: u vs v   
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (6,3)
# fig4, ax = plt.subplots(nrows=1, ncols=1)
 
# ax.set_xlabel('u  [ pA ]'); ax.set_ylabel('v  [ mV ]')

# ax.set_title('C = %0.0f' %C+ '    k = %0.2f' %k +
#              '    v$_r$ = %0.0f' %vr + '    v$_t$ = %0.0f' %vt + 
#              '\n a = %0.2f' %a + '   b = %0.1f' %b +
#             '   c = %0.0f' %c + '    d = %0.0f' %d +
#             '   I$_0$ = %0.1f' %I0 , fontsize = 11)

# ax.grid()

# ax.plot(u,v,'b', lw = 1)
# ax.plot(u[0],v[0],'bo', ms = 6)

# fig4.tight_layout()
# fig4.savefig('a4.png')


# #%%   Response to a ramp input: ISI 
# peaks = find_peaks(v)
# Peaks = peaks[0]
# L = len(Peaks)
# tPeaks = t[Peaks]
# tPeaks2 = tPeaks[1:L]
# tPeaks1 = tPeaks[0:L-1]
# ISI = tPeaks2 - tPeaks1
# f = 1000/ISI
# tPEAKS = (tPeaks2+tPeaks1)/2
# index = (Peaks[0:L-1] + Peaks[1:L])/2
# index = index.astype(int)
# IPEAKS = Iext[index]

# print('I0 = %0.1f' %I0)
# q = np.mean(f); print('mean f = %0.2f' %q)
# q = np.mean(ISI); print('mean ISI = %0.0f' %q)

# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (6,3)
# fig5, ax = plt.subplots(nrows=1, ncols=1)
 
# ax.set_xlabel('$I_{ext}$  [ pA ]'); ax.set_ylabel('f  [ Hz ]')

# ax.set_title('C = %0.0f' %C+ '    k = %0.2f' %k +
#              '    v$_r$ = %0.0f' %vr + '    v$_t$ = %0.0f' %vt + 
#              '\n a = %0.2f' %a + '   b = %0.1f' %b +
#             '   c = %0.0f' %c + '    d = %0.0f' %d 
#             , fontsize = 11)

# ax.grid()
# ax.set_xlim([0,max(IPEAKS)])
# ax.set_ylim([0,1.1*max(f)])

# ax.plot(IPEAKS,f,'go', ms = 3)

# fig5.tight_layout()
# fig5.savefig('a5.png')
         
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


