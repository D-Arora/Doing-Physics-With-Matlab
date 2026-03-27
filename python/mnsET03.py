# -*- coding: utf-8 -*-
'''

mnsET02.py.py      MARCH 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

MORRIS - LECAR MODEL OF A NEURON

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mnsET03.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, log, linspace, zeros, array, ones, sqrt, exp 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.linalg import eig
import time
import matplotlib.animation as animation

tStart = time.time()
plt.close('all')


#%% FUNCTIONS
def VDOT(t,v,u,Iext):
    Vdot  = (ICa + IK + IL + I)/C
    return Vdot    
def UDOT(t,v,u,Iext):
    Udot  = phi*(UI - u)/tau
    return Udot   

def RK(t,v,u,I):
    k1 = VDOT(t,v,u,I)
    g1 = UDOT(t,v,u,I)
    
    k2 = VDOT(t+0.5*h,v+k1*h/2,u+g1*h/2,I)
    g2 = UDOT(t+0.5*h,v+k1*h/2,u+g1*h/2,I)
    
    k3 = VDOT(t+0.5*h,v+k2*h/2,u+g2*h/2,I)
    g3 = UDOT(t+0.5*h,v+k2*h/2,u+g2*h/2,I)
    
    k4 = VDOT(t+1.0*h,v+k3*h,u+g3*h,I)
    g4 = UDOT(t+1.0*h,v+k3*h,u+g3*h,I)
    
    vNew = v + h*(k1 + 2*k2 + 2*k3 + k4)/6
    uNew = u + h*(g1 + 2*g2 + 2*g3 + g4)/6
    
    return [vNew,uNew]

#%%
# Hopf data
C = 20.0       # Membrane capacitance
gCa = 4.4; ECa = 120.0
gK = 8.0;  EK = -84.0
gL = 2.0;  EL = -60.0
v1 = -1.2; v2 = 18.0; v3 = 2.0; v4 = 30.0
phi = 0.04 


# Saddle Node on a LimitCycle (SNLC) data
#phi = 0.067
#gCa = 4
#v3 = 12
#v4 = 17.4

# Saddle Node Homoclinic data
phi = 0.23
v4 = 17.4

I = 100
col = [0,0,1]
tMax = 200

nT = 9999

t = linspace(0,tMax,nT)
h = t[2] - t[1]
v = zeros(nT); u = zeros(nT)
mI = zeros(nT); uI = zeros(nT)
iCa = zeros(nT); iK = zeros(nT); iL = zeros(nT)

# Initial Conditions
#v[0],u[0] = 0.2,0.18
#v[0],u[0] = -9,0.08
#v[0],u[0] = -60,0.05
#v[0], u[0] = 0,0
v[0], u[0] = 0, 0.2

#%% Runge-Kutta: solve ODEs
for n in range(nT-1):
    V = v[n]; U = u[n]
    mI[n] = 0.5 * (1 + np.tanh((V - v1) / v2))
    uI[n] = 0.5 * (1 + np.tanh((V - v3) / v4)); UI = uI[n]
    tau = 1 / ( np.cosh((V - v3) / (2 * v4)) )
    ICa = -gCa * mI[n] * (V - ECa)
    IK = -gK * U * (V - EK)
    IL = -gL * (V - EL)
    iCa[n] = ICa; iK[n]= IK; iL[n] = IL
        
    v[n+1], u[n+1] = RK(t[n],v[n],u[n],I)

mI[-1] = mI[-2]; uI[-1] = uI[-2]

#%%  NULLCLINES
vN = linspace(-70,20,599)
# u nullcline
uNu = 0.5 * (1 + np.tanh((vN - v3) / v4))
MI = 0.5 * (1 + np.tanh((vN - v1) / v2))
# v nullcline
uNv = (I - gL*(vN - EL) - gCa*MI*(vN - ECa) ) / (gK*(vN - EK))
 
# Find indices where sign changes --> steady values for v and u
vu = uNu - uNv  

index = np.where(np.diff(np.sign(vu)))[0]
vSS = vN[index]
uSS = 0.5 * (1 + np.tanh((vSS - v3) / v4))

print(f"Indices of zero crossings: {index}")
print("Values near crossings")
print('vSS')
print(np.round(vSS,2))
print('uSS')
print(np.round(uSS,3))      

#xxx
#%%   FIG 1:   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('v [ mV ]')
ax.grid()
ax.set_title('Membrane potential    I$_{ext}$ = %0.1f pA' %I, fontsize = 12)

ax.plot(t,v,lw = 2, color = col)

ax.plot([0,tMax],[vSS,vSS],'k',lw = 1) 

#ax.set_xlim([0,50])
fig1.tight_layout()
fig1.savefig('a1.png')

#%%   FIG 2:   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_title('Recovery variable    I$_{ext}$ = %0.1f' %I, fontsize = 12) 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('u')
ax.grid()
ax.plot(t,u,'k',lw = 2) 
fig2.tight_layout()
fig2.savefig('a2.png')

#%%   FIG 3:   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)
ax.set_title('phase portrait    I$_{ext}$ = %0.1f pA' %I, fontsize = 12) 
ax.set_xlabel('v  [ mV ]'); ax.set_ylabel('u')
ax.grid()
ax.plot(v,u,'b',lw = 2) 
ax.plot(vN,uNu,'m',lw = 1, label = 'u$_{null}$')
ax.plot(vN,uNv,'r',lw = 1, label = 'v$_{null}$')
ax.plot(v[0],u[0], 'go', ms = 6)
ax.legend(fontsize = 10)
fig3.tight_layout()
fig3.savefig('a3.png')

#%%   FIG 4:   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_title('Gate variables    I$_{ext}$ = %0.1f' %I, fontsize = 12) 
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('mI, uI')
ax.grid()
ax.plot(t,mI,'r',lw = 2, label = 'm$_I$')
ax.plot(t,uI,'b',lw = 2, label = 'u$_I$') 
ax.legend(fontsize = 10)
fig4.tight_layout()
fig4.savefig('a4.png')

#%%   FIG 5:   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('currents [ pA]')
ax.set_title('Ion Currents    I$_{ext}$ = %0.1f  pA' %I, fontsize = 12)
ax.grid()
ax.plot(t,iCa,'r',lw = 2, label = 'i$_{Ca}$')
ax.plot(t,iK,'b',lw = 2,label = 'i$_K$') 
ax.plot(t,iL,'m', lw =2,label = 'i$_L$')
ax.legend(fontsize = 10)
fig5.tight_layout()
fig5.savefig('a5.png')


#%% Fig 6   Vector field
X = linspace(-70,100,12); Y = linspace(0,1,12)
xx,yy = np.meshgrid(X,Y)
tuI = 0.5 * (1 + np.tanh((xx - v3) / v4))
tmI = 0.5 * (1 + np.tanh((xx - v1) / v2))
ttau = 1 / ( np.cosh((xx - v3) / (2 * v4)) )

xxDot = -gCa*tmI*(xx-ECa)-gK*yy*(xx-EK)-gL*(xx-EL)+I
yyDot = (phi/ttau)*(tuI-yy)

XXDot = xxDot/(sqrt(xxDot**2 + yyDot**2))
YYDot = yyDot/(sqrt(xxDot**2 + yyDot**2)) 

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,4)
fig6, ax = plt.subplots(nrows=1, ncols=1)
ax.set_title('Vector field    I$_{ext}$ = %0.1f pA' %I, fontsize = 12) 
ax.set_xlabel('v  [ mV ]'); ax.set_ylabel('u')
#ax.quiver(xx,yy,XXDot,YYDot)
ax.streamplot(xx,yy,xxDot,yyDot, density = [0.7,0.7])
ax.plot(vN,uNu,'m',lw = 1, label = 'u$_{null}$')
ax.plot(vN,uNv,'r',lw = 1, label = 'v$_{null}$')

ax.plot(v[0],u[0], 'go', ms = 6)
ax.plot(v,u,'b',lw = 2)

ax.plot(vSS,uSS, 'ro', ms = 7) 
ax.legend(fontsize = 10)
ax.set_ylim([-0.05,0.6])
fig6.tight_layout()
fig6.savefig('a6.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


