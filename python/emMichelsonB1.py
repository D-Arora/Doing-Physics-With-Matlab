# -*- coding: utf-8 -*-
"""
emMichelsonB1.py          mar 2025

COMPUTATIONAL OPTICS
MICHELSON INTERFEROMETER: TWO POINT SOURCES
 
Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/emMichelson..pdf
    

"""

#%%
import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array 
from scipy import pi, sqrt
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation, PillowWriter 

#%%
N = 999        # Grid points

wL = 700e-9; col = [1,0,0]     # red wavelength [m]
#wL = 550e-9; col = [0,1,0]    # green wavelength [m]
#wL = 480e-9; col = [0,0,1]    # blue wavelength [m]


k = 2*pi/wL    # Propagation constant  [1/m]

# Detector space - observation screen
L = 700e-9
zD = 60*L
XD = 80*L
xD = linspace(-XD,XD,N)

# Two point sources: Source 1 at Origin (0,0,0) 
x1 = 0; x2 = 0
z1 = 0; z2 = -10.0*L

# Calculation of electric fields and screen intensity S  [a.u.]
r1 = ( (xD - x1)**2 + (zD-z1)**2 )**0.5
r2 = ( (xD - x2)**2 + (zD-z2)**2 )**0.5
E1 = exp(1j*k*r1)/r1
E2 = exp(1j*(k*r2+pi))/r2
E = E1+E2
S = np.real(np.conj(E)*E)


#%%   GRAPHICS
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)


xP = xD/wL; yP = S = S/max(S)
ax.plot(xP,yP,color = col,lw = 2)
ax.grid()
ax.set_ylabel('S  [a.u.]')
ax.set_xlabel('x$_D / \lambda$')

s1 = wL*1e9; s2 = x2/wL; s3 = z2/wL
ax.set_title('$\lambda = $ %0.0f nm' %s1
             +  '  $z_D$ =  %0.2e  ' %zD
             + '   $x_2 / \lambda = $ %0.1f'  %s2 
             + '   $z_2 / \lambda = $ %0.1f'  %s3 )


fig1.tight_layout()
fig1.savefig('a1.png') 


#%%   ANIMATION

z  = linspace(-10*L,-1*L,100)
z = np.flip(z)
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,4)
fig3, ax = plt.subplots(1)

ax.set_xlim([-80,80])
ax.set_ylim([0,1.2])
ax.set_xlabel('x$_D / \lambda$')
txt = ax.text(-75,1.0, 'x',fontsize=12)
line1, = ax.plot([], [], lw = 2, color = col)     
plt.yticks([])   
fig3.tight_layout()


def init():  
    line1.set_data([], [])
    txt.set_text('')
    return  line1, txt, 
   
def animate(n):
       r2 = ( (xD - x2)**2 + (zD-z[n])**2 )**0.5
       E2 = exp(1j*(k*r2+pi))/r2
       E = E1+E2
       S = np.real(np.conj(E)*E)
       xP = xD/wL; yP = S = S/max(S)
       u1 = xP; v1 = yP
       line1.set_data([u1], [v1]) 
       time.sleep(0.2)
       s = z[n]/wL; s1 = abs(s); txt.set_text('$z_2 / \lambda$ = %2.2f   ' % s + ' N = %2.0f' %s1  )  
       return   line1, txt,

anim = FuncAnimation(fig3, animate, init_func = init, 
                      frames = 100, interval = 2, blit = True, repeat = False)

#  anim.save('agMIB1.gif', fps = 5)




