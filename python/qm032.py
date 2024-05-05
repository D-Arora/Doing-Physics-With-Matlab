# -*- coding: utf-8 -*-
"""
# -*- coding: utf-8 -*-

qm007.py    April 2024

QUANTUM MECHANICS
Finite Difference Time Development Method: Animation
     [1D] Schrodinger Equation:
     Free particle: wavepacket propagation in a parabolic potential well

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm007.htm


"""

# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, real, imag
import cmath
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation, PillowWriter 

#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [ J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
me = 9.10938356e-31            # Electron mass [kg]

xMin = -0.6e-9; xMax = 0.6e-9     # X domain  [m]
E0 = 80                   # beam particle energy  [eV]
U0 = 50                  # height of potential step [eV]
N = 599                           # Grid points

            
#%% COMPUTATIONS
x = np.linspace(xMin,xMax,N)
dx = x[2] - x[1]
E = e*E0                          # beam particle energy  [J]                                 
U = e*U0                          # height of potential step [J]
w = E/hbar                        # angular frequency [rad/s]
T = 2*pi/w                        # oscillation period [s]
k1 = sqrt(2*me*E)/hbar            # propagation constant x < 0    [1/m]
k2 = cmath.sqrt(2*me*(E-U))/hbar       # propagation constant x > 0    [1/m]
L1 = 2*pi/k1                      # wavelength x < 0    [m]
L2 = 2*pi/k2                      # wavelength x > 0    [m]


#%%
# WAVEFUNCTIONS:  region 1 and region 2
A = 1
B = (k1 - k2) / (k1 + k2)
C = 2*k1 / (k1 + k2)

x1 = np.linspace(xMin,0,N)
y1 = A*exp(1j*k1*x1)
y1R = real(y1)
y1I = imag(y1)

z1 = B*exp(-1j*k1*x1)
z1R = real(z1)
z1I = imag(z1)

psi1 = y1 + z1
psi1R = real(psi1)
psi1I = imag(psi1)

x2 = np.linspace(0,xMax,N)
y2 = C*exp(1j*k2*x2)
psi2R = real(y2)
psi2I = imag(y2)

nT = 120
P = 2*pi/w
t = linspace(0,5*T,nT)
T = exp(-1j*w*t)


#%% GRAPHICS
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,4)
fig, ax = plt.subplots(1)
fig.subplots_adjust(top = 0.92, bottom = 0.20, left = 0.20,\
                    right = 0.92, hspace = 0.20,wspace=0.2)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('$\psi_R$',color= 'black',fontsize = 16)
ax.set_xlabel('x  [nm]',color = 'black')
ax.set_xlim([-0.6, 0.6])
ax.set_ylim([-1.20, 1.2])
#ax.set_xticks(np.arange(0,101,20))
#ax.set_yticks(np.arange(-20,81,20))
#txt.set_text('time = %.2f' % t[n] + ' s')  
ax.set_title('E0 = %2.0f' % E0 + '    U0 = %2.0f' % U0, fontsize = 12)

line1, = ax.plot([], [], 'b', lw = 2)     

line2, = ax.plot([], [], 'b', lw = 2)

#ax.grid('visible')

def init():  
   line1.set_data([], [])
   line2.set_data([], [])
   return  line1, line2
   
def animate(n):
     u1 = x1*1e9
     v1 = real(psi1R*T[n]/max(psi1R))
     line1.set_data([u1], [v1]) 
     
     u2 = x2*1e9
     v2 = real(y2*T[n]/max(psi1R))
     if E0 < U0:
        v2 = real(psi2R*T[n]/max(psi1R))
     line2.set_data([u2], [v2]) 
     
     time.sleep(0.1)
     return   line1, line2

anim = FuncAnimation(fig, animate, init_func = init, 
                     frames = nT, interval = 1, blit = True, repeat = False)

anim.save('agA.gif', fps = 5)
#anim.save('agB.gif', fps = 5)

