# -*- coding: utf-8 -*-

"""
op001.py    oct 2024

COMPUTATIONAL OPTICS

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/op001.htm


INTERFERENCE

TRAVELLING MONOCHROMATIC WAVE

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, real

tStart = time.time()


#%%
# Input wavelength  [m]   >>>
wL = 500e-9

# Model parameters
c = 3e8        # speed of light in vaccum  [m/s]
k = 2*pi/wL    # propagation constant  [rad/m] 
f = c/wL       # frequency  [Hz]  
w = 2*pi*f     # angular frequency  [rad/s]
T = 1/f        # period  [s]
A = 1          # amplitude
Nt = 99        # number of steps 
Nz = 999
tMax = 3*T; zMax = 5*wL
t = linspace(0,tMax,Nt)
z = linspace(0,zMax,Nz)
 
tG, zG = np.meshgrid(t,z)

# Wavefunction: --> (-)   <-- (+)
u = A*exp(1j*(k*zG + w*tG))

# Create a figure and axes
fig = plt.figure(figsize=(5,3))
fig.subplots_adjust(top=0.85, bottom = 0.20, left = 0.15,\
                    right = 0.95, hspace = 0.60,wspace=0.5)
ax1 = plt.subplot(1,1,1)   
ax1.set_xlim(( 0, zMax*1e9))            
ax1.set_ylim((-1.1*A, 1.1*A))
ax1.set_xlabel('z  [nm]')
ax1.set_ylabel('u  [a.u.]')
txt_title = ax1.set_title('')
line1, = ax1.plot([], [], 'blue', lw=2)     
plt.grid('visible')
#fig.tight_layout()

# Animation
def drawframe(n):
    xf = z*1e9
    yf = np.real(u[:,n])
    line1.set_data(xf, yf)
    s = t[n]/T
    txt_title.set_text('t/T = {:.2f}'.format(s))
    
    ax1.set_xlim(( 0, zMax*1e9))            
    ax1.set_ylim((-1.1*A, 1.1*A))
    return line1

anim = animation.FuncAnimation(fig, drawframe, frames=Nt, interval=100, blit=False, \
                               repeat = False)


#%% Position and time plots
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=2)
#fig1.subplots_adjust(top = 0.94, bottom = 0.15, left = 0.180,\
#                    right = 0.92, hspace = 0.36,wspace=0.40)
R = 0
ax[R].set_xlabel('z  [nm]',fontsize = 12)
ax[R].set_ylabel('u  [a.u.]',fontsize = 12)
ax[R].xaxis.grid()
ax[R].yaxis.grid()
xP = z*1e9; yP = real(u[:,0])
ax[R].plot(xP,yP,'b',lw = 2)

R = 1
ax[R].set_xlabel('t  [fs]',fontsize = 12)
ax[R].set_ylabel('u  [a.u.]',fontsize = 12)
ax[R].xaxis.grid()
ax[R].yaxis.grid()
xP = t*1e15; yP = u[0,:]
ax[R].plot(xP,yP,'r',lw = 2)
fig1.tight_layout()  

  
#%% Save figures    
# anim.save("ag_A.gif", dpi=250, fps=10)
# fig1.savefig('a1.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.3f s' %tExe)