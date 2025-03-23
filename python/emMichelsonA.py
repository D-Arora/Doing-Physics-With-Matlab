# -*- coding: utf-8 -*-
"""
emMichelsonA.py          mar 2025

COMPUTATIONAL OPTICS: MICHELSON INTERFEROMETER
   source: plane waves
   Plot detector screen intensity
   Animation of screen intensity as separation of mirror M2 chnages

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/emMichelson.pdf
"""

#%%
import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array 
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation, PillowWriter 


#%%   Calculate screen intensity
num = 999     # Grid points
Lmax = 2      # Detector screen dimensions
L = linspace(-Lmax,Lmax, num)
SD = (sin(pi*L))**2     # Detector screen intensity for plane wave sources


#%%  [1D] plot intensity
plt.rcParams["figure.figsize"] = (4,2.5)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('2 $\Delta$d / $\lambda$ ',fontsize = 12)
ax.set_ylabel('$S_D$  [ a.u ] ',fontsize = 12)
ax.grid()
xP = L; yP = SD
ax.plot(xP,yP,'r',lw = 2)
fig1.tight_layout()


#%%  [2D] plot intensity
plt.rcParams["figure.figsize"] = (6,7)
fig2, ax = plt.subplots(nrows=3, ncols=3)

def circle(R,C,col):
    ax[R,C].axis('off')
    ax[R,C].plot(0,0,'o',ms = 80, color = [col,0,0])
    ax[R,C].set_title('2 $\Delta$d / $\lambda$ = %0.3f' %col)
    fig2.tight_layout()
    
R = 0; C = 0; col = 0/8; circle(R,C,col)
R = 0; C = 1; col = 1/8; circle(R,C,col)
R = 0; C = 2; col = 2/8; circle(R,C,col)

R = 1; C = 0; col = 3/8; circle(R,C,col)
R = 1; C = 1; col = 4/8; circle(R,C,col)
R = 1; C = 2; col = 5/8; circle(R,C,col)

R = 2; C = 0; col = 6/8; circle(R,C,col)
R = 2; C = 1; col = 7/8; circle(R,C,col)
R = 2; C = 2; col = 8/8; circle(R,C,col)


#%% ANIMATION   intensity

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2)
fig3, ax = plt.subplots(1)

ax.set_xlim([-1,1])
ax.set_ylim([-1,1])

ax.axis('off')
n = 0; col = [0,0,0]
txt = ax.text(-1,0.5, 'x',fontsize=12)
line1, = ax.plot([], [], 'o',ms = 80, color = col)     
    
fig3.tight_layout()


def init():  
   line1.set_data([], [])
   txt.set_text('')
   return  line1, txt, 
   
def animate(n):
     L = (2/(199))*n
     col = (sin(pi*L))**2
     line1, = ax.plot(0,0,'o',ms = 80, color = [col,0,0]) 
     time.sleep(0.1)
     txt.set_text('$\Delta$ = %2.2f' % L )  
     return   line1, txt,

anim = FuncAnimation(fig3, animate, init_func = init, 
                     frames = 200, interval = 10, blit = True, repeat = False)

# anim.save('agA.gif', fps = 5)


#%% Plant Growth Calculation
wL = 632.8e-9;        # wavelength [m]
dt = 8.3              # Time interval  [h]
nf = 3420;            # Number of fringes
dP = (nf * wL /2) * 1e3;   # distance moved by plant on mirror 2 [mm]
dPdt = dP/dt;         # rate of growth  [mm/h]
  
print('Inputs  ')
print('   wavelength = %3.1e  m ' % wL)
print('   time interval = %3.1f  h' % dt)
print('   fringes  = %3.0f ' %nf)
print('Outputs  ')
print('   growth distance = %3.5f  mm' % dP)
print('   rate of growth  = %3.5f  mm/h ' % dPdt)


#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')