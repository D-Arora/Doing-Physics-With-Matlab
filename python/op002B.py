# -*- coding: utf-8 -*-

"""
op002B.py    oct 2024

COMPUTATIONAL OPTICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/op002.htm


INTERFERENCE: BEATS



"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag, real
import random
tStart = time.time()


#%%
# Inputs: wavelength [m], amplitude, initial phase angle [rad]   >>>
wL = zeros(3); wL[0] = 500e-9; wL[1] = 560e-9; wL[2] = 300e-9
A = zeros(3); A[0] = 1; A[1] = 1; A[2] = 0
phi = zeros(3); phi[0] = 0; phi[1] = 0; phi[2] = 0
Amax = 2.2      # max amplitude for plots
Nt = 299        # number of steps 
Nz = 999        # Z axis grid
tEnd = 15       # number of periods for simulation for time grid
zEnd = 20       # number of wavelengths for z grid

#%% SETUP
# Model parameters
c = 3e8        # speed of light in vaccum  [m/s]
k = 2*pi/wL    # propagation constant  [rad/m] 
f = c/wL       # frequency  [Hz]  
w = 2*pi*f     # angular frequency  [rad/s]
T = 1/f        # period  [s]
tMax = tEnd*T[0]; zMax = zEnd*wL[0]
t = linspace(0,tMax,Nt)
z = linspace(0,zMax,Nz)
 
tG, zG = np.meshgrid(t,z)

# Wavefunctions: --> (-)   <-- (+)
u1 = A[0]*exp(1j*(k[0]*zG - w[0]*tG  + phi[0]))
u2 = A[1]*exp(1j*(k[1]*zG - w[1]*tG  + phi[1]))
u3 = A[2]*exp(1j*(k[2]*zG - w[2]*tG  + phi[2]))

u = u1 + u2 + u3

#%% Console output
wLs = 1e9*wL; Ts = 1e15*T
print('wL [nm]   A     f [Hz]     T [fs]')
for m in range(3):
    print(' %2.0f ' %wLs[m] + '    %2.0f' %A[m] + '   %2.3e' %f[m]  
          + '   %  2.3f' %Ts[m] ) 

fBeat = abs(f[0]-f[1])
TBeat = 1e15/fBeat

fFast = (f[0]+f[1])/2
TFast = 1e15/fFast

print('  ')
print('f_beats [Hz]   T_beats [fs]   f_fast [Hz]   T_fast [fs]')
print('%2.3e' %fBeat + '      %2.3f' %TBeat + 
      '         %2.3e' %fFast + '     %2.3f' %TFast)
    
#%%
# Create a figure and axes
fig = plt.figure(figsize=(5,3))
fig.subplots_adjust(top=0.85, bottom = 0.20, left = 0.15,\
                    right = 0.95, hspace = 0.60,wspace=0.5)
ax1 = plt.subplot(1,1,1)   
ax1.set_xlim(( 0, zMax*1e9))            
ax1.set_ylim((-Amax, Amax))
ax1.set_xlabel('z  [nm]')
ax1.set_ylabel('u  [a.u.]')
line1, = ax1.plot([], [], 'blue', lw=2)  
line2, = ax1.plot([], [], 'k', lw=1) 
line3, = ax1.plot([], [], 'k', lw=1) 
line4, = ax1.plot([], [], 'k', lw=1)     
plt.grid('visible')
fig.tight_layout()

# Animation
def drawframe(n):
    xf = z*1e9
    yf = np.real(u[:,n])
    line1.set_data(xf, yf)
    yf = np.real(u1[:,n])
    line2.set_data(xf, yf)
    yf = np.real(u2[:,n])
    line3.set_data(xf, yf)
    yf = np.real(u3[:,n])
    line4.set_data(xf, yf)
    ax1.set_xlim(( 0, zMax*1e9))            
    ax1.set_ylim((-Amax, Amax))
    return line1,line2,line3,line4,

anim = animation.FuncAnimation(fig, drawframe, frames=Nt, interval=50
                               
                               , blit=False, \
                               repeat = False)


#%% Position and time plots

plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
#fig1.subplots_adjust(top = 0.94, bottom = 0.15, left = 0.180,\
#                    right = 0.92, hspace = 0.36,wspace=0.40)

ax.set_xlabel('z  [nm]',fontsize = 12)
ax.set_ylabel('u  [a.u.]',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
xP = z*1e9; yP = real(u[:,0])
ax.plot(xP,yP,'b',lw = 2)
yP = real(u1[:,0])
ax.plot(xP,yP,'k',lw = 1)
yP = real(u2[:,0])
ax.plot(xP,yP,'k',lw = 1)
yP = real(u3[:,0])
ax.plot(xP,yP,'k',lw = 1)
fig1.tight_layout()

#%%
plt.rcParams["figure.figsize"] = (6,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
#fig1.subplots_adjust(top = 0.94, bottom = 0.15, left = 0.180,\
#                    right = 0.92, hspace = 0.36,wspace=0.40)

ax.set_xlabel('t  [fs]',fontsize = 12)
ax.set_ylabel('u  [a.u.]',fontsize = 12)
ax.xaxis.grid()
ax.yaxis.grid()
xP = t*1e15; yP = np.real(u[0,:])
ax.plot(xP,yP,'r',lw = 2)
yP = np.real(u1[0,:])
ax.plot(xP,yP,'k',lw = 1)
yP = np.real(u2[0,:])
ax.plot(xP,yP,'k',lw = 1)
yP = np.real(u3[0,:])
ax.plot(xP,yP,'k',lw = 1)
fig2.tight_layout()  




#%% RANDOM AND COHERENT SOURCES
nR = 100
phiR = zeros(nR)
uR = zeros(nR) + 1j*zeros(nR)
for m in range(nR):
    phiR[m] = 2*pi*random.random() 
    uR[m] = exp(1j*phiR[m])

IR = np.conj(uR)*uR
Itot = sum(IR)
#print(Itot)
  
#%% Save figures    
# anim.save("ag_A.gif", dpi=250, fps=8)
# fig1.savefig('a1.png')
# fig3.savefig('a2.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)