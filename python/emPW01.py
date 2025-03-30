# -*- coding: utf-8 -*-
"""
emPW01.py          mar 2025

COMPUTATIONAL OPTICS: SUPERPOSITION OF PLANE WAVES
     pLANES WAVES TRAVEL IN THE +z DIRECTION

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/emPW01.htm
"""

#%%
import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array, real, imag, conj 
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation, PillowWriter 


#%%   INPUTS
num = 999     # Grid points
N = 2; wL = zeros(N); n = zeros(N)
# Wavelengths wL [m];  wL[0] < wL[1]
wL[0] = 20e-3;  wL[1] = 25e-3  
# refractive indices  n:  n[0] > n[1]
# longer wavelength must have a lower value for refracticve its refractive index 
n[0]  = 1.1;  n[1] = 1.100     


#%%   CALCULATIONS
c = 3e8           # speed of light
k = 2*pi/wL        # propagation constant
w = c*k/n          # omega: angular frequency
T = 2*pi/w[0]      # period wave 1
vP = w/k            # phase velocities
vG = (w[1] - w[0]) / (k[1]- k[0])

L = 20*wL[0]    # Z axis range
z = linspace(0,L,num)
F = 200; t = linspace(0,20*T,F)


# Spatial electric field at time ts
ts = 0*T      # ts = 0 for animation
E0z = exp(1j*k[0]*z)*exp(-1j*w[0]*ts)
E1z = exp(1j*k[1]*z)*exp(-1j*w[1]*ts)
Ez = E0z + E1z

S = 0*Ez        # intensity

#%%
wp = (w[0]+w[1])/2; kp = (k[0]+k[1])/2
wg = (w[0]-w[1])/2; kg = (k[0]-k[1])/2
vG = wg/kg
vR = wp/kp

#%%  [1D]    PLOT ELECTRIC FIELS AT TIME t = 0 
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('z  [m]',fontsize = 12)
ax.set_ylabel('E(0) [ a.u ] ',fontsize = 12)
q = 1e3*wL; ax.set_title('t = %0.2e' %ts+
                         '   $\lambda_0$ = %0.0f mm   ' %q[0] + 
                         '$n_0$ = %0.2f' %n[0] +
                         '   $\lambda_1$ = %0.0f mm' %q[1] +
                         '   $n_1$ = %0.2f   ' %n[1]
                         , fontsize = 10)
ax.grid()
xP = z
yP = E0z; ax.plot(xP,yP,'r',lw = 1,label = 'E$_{0z}$')
yP = E1z; ax.plot(xP,yP,'b',lw = 1,label = 'E$_{1z}$') 
yP = Ez;  ax.plot(xP,yP,'k',lw = 2,label = 'E$_z$') 

#ax.legend()
fig1.tight_layout()


#%%  CONSOLE OUTPUT
print('wL0 = %0.3f m' %wL[0] + '   wL1 =  %0.3f m' % wL[1])
print('n0 = %0.3f' %n[0] + '   n1 =  %0.3f m' % n[1])
print('vP0 = %0.3e m/s' %vP[0] + '   vP1 =  %0.3e m/s' %vP[1])
print('vG = %0.3e m/s' %vG)
print('vR = %0.3e m/s' %vR)

#%% ANIMATION   travelling plane waves

E0 = E0z; E1 = E1z

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,5)
fig3, ax = plt.subplots(1)

ax.set_xlim([0,L])
ax.set_ylim([-0.20,11.5])

ax.axis('off')

line1, = ax.plot([], [], 'r',lw = 2)     
line2, = ax.plot([], [], 'b',lw = 1)  
line3, = ax.plot([], [], 'k',lw = 2)  
line4, = ax.plot([], [], 'm',lw = 2) 
line5, = ax.plot([], [], 'mo',ms = 8) 
line6, = ax.plot([], [], 'ro',ms = 8)  
line7, = ax.plot([], [], 'bo',ms = 8)  
line8, = ax.plot([], [], 'm',lw = 1)      
fig3.tight_layout()


def init():  
    return  line1, line2, line3, line4, line5, line6, line7,line8,
 
   
def animate(q):
       E0 = E0z*exp(-1j*w[0]*t[q])
       u1 = z; v1 = real(E0)+10
       line1.set_data([u1], [v1]) 
       
       E1 = E1z*exp(-1j*w[1]*t[q])
       u1 = z; v1 = real(E1)+10
       line2.set_data([u1], [v1]) 
       
       u1 = z
       v1 = real(E0+E1)+6.5
       line3.set_data([u1], [v1])
       
       E = E0 + E1
       S = conj(E)*E
       v1 = S
       line4.set_data([u1], [v1])
       
       v1 = real(E)**2
       line8.set_data([u1], [v1])
       
       u1 = vG*t[q]; v1 = 4
       line5.set_data([u1], [v1])
       
       u1 = vP[0]*t[q]; v1 = 11
       line6.set_data([u1], [v1])
       
       u1 = vP[1]*t[q]; v1 = 11
       line7.set_data([u1], [v1])
           
       time.sleep(0.05)
       return   line1, line2, line3, line4, line5, line6, line7, line8

anim = FuncAnimation(fig3, animate, init_func = init, 
                     frames = F, interval = 2, blit = True, repeat = False)

# anim.save('agPW01.gif', fps = 10)


#%%
# fig1.savefig('a1.png')
