# -*- coding: utf-8 -*-
"""
qm042.py    May 2024

QUANTUM MECHANICS
    Vibrations of a guitar string
    Solving the time independent wave equation for standing waves on a string
    by finding the eigenvalues and eigenfunctions
    Boundar conditons: fixed - free

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation 
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm042an.htm

    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm042.pdf


"""

#%%  IMPORTS
import numpy as np
from scipy.sparse import diags #Allows us to construct our matrices
from scipy.sparse.linalg import eigsh #Solves the Eigenvalue problem
import matplotlib.pyplot as plt
from scipy.linalg import *
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, real, imag
import time
from matplotlib.animation import FuncAnimation, PillowWriter 

tStart = time.time()


#%% INPUTS
N = 99          # grid points
L = 1           # length of guitar string  0.640  [m]
F = 100         # string tension  [N]
mu = 5.0e-3     # linear density of string [kg/m]
M = 6           # Number of modes returned from eigsh 
nM = 4          # Mode number for animation

#%% COMPUTATIONS
x = linspace(0,L,N)      # x domain  [m]
dx = x[2] - x[1]          # dx [m]        

#  Second Derivative matrix: default fixed ends at x = 0 and x = L
off = np.ones(N-1)
A = -2*np.eye(N) + np.diag(off,1) + np.diag(off,-1) 

A[-1,-1] = -1     # free end at x = L  
#A[0,0]   = -1     # free end at x = 0

# Eigenvalues and eigenfunctions (eigenvectors)
ev, ef = eigsh(A, which="SM", k = M)
ef = ef / np.amax(ef)

#%% Normal modes
wL = np.sqrt(-4*np.pi**2*dx**2/ev)   # wavelengths  [m]
v = sqrt(F/mu)                       # velocity  [m/s]
f = v/wL                             # frequency  [Hz]
T = 1/f                              # period  [s]    
w = 2*pi*f                           # angular frequency  [rad/s]

#%% Console display
print(' ')
print('L = %2.1f m' % L + '  F = %2.1f N' % F + '   mu = %2.1e kg/m' %mu)
print('v = %2.2f m/s' % v)

print('  n    lambda [m]    f  [Hz]     fn/f1')
fR = zeros(M)
for c in range(M):
    sc = 2*c+1; s = M-1-c
    fR[s] = f[s]/f[M-1]
    print(' %2.0f' % sc + '       %2.2f' % wL[s]
              + '        %2.2f'    % f[s] + '      %2.2f'    % fR[s]  )


#%% GRAPHICS

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7)
fig1, axes = plt.subplots(nrows=3, ncols=2)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
                    right = 0.95, hspace = 0.36,wspace=0.40)

plt.suptitle('Normal modes')

xP = x

R = 0; C = 0; mode = 5   
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_1$ = %2.2f m' % wL[mode], fontsize = 10)
#yP = zeros(N+2); yP[1:-1] = ef[:,mode]
#xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
yP = ef[:,mode]
axes[R,C].plot(xP,yP, 'blue')

R = 0; C = 1; mode = 4   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_3$ = %2.2f m' % wL[mode], fontsize = 10)
#yP = zeros(N+2); yP[1:-1] = ef[:,mode]
#xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
yP = ef[:,mode]
#axes[R,C].plot(xP,yP, 'blue')
axes[R,C].plot(xP,yP, 'blue')

R = 1; C = 0; mode = 3   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_5$ = %2.2f m' % wL[mode], fontsize = 10)
#yP = zeros(N+2); yP[1:-1] = ef[:,mode]
#xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x; 
yP = ef[:,mode]
axes[R,C].plot(xP,yP, 'blue') 

R = 1; C = 1; mode = 2   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_7$ = %2.2f m' % wL[mode], fontsize = 10)
#yP = zeros(N+2); yP[1:-1] = ef[:,mode]
#xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
yP = ef[:,mode]
axes[R,C].plot(xP,yP, 'blue')

R = 2; C = 0; mode = 1   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_9$ = %2.2f m' % wL[mode], fontsize = 10)
#yP = zeros(N+2); yP[1:-1] = ef[:,mode]
#xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
yP = ef[:,mode]
axes[R,C].plot(xP,yP, 'blue')

R = 2; C = 1; mode = 0   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_{11}$ = %2.2f m' % wL[mode], fontsize = 10)
#yP = zeros(N+2); yP[1:-1] = ef[:,mode]
#xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
yP = ef[:,mode]
axes[R,C].plot(xP,yP, 'blue')

fig1.savefig('a1.png')


#%% STANDING WAVE VIEW
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (7,7)
fig1, axes = plt.subplots(nrows=3, ncols=2)
fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
                    right = 0.95, hspace = 0.36,wspace=0.40)

plt.suptitle('Normal modes')

R = 0; C = 0; mode = 5   
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_1$ = %2.2f m' % wL[mode], fontsize = 10)
yP = zeros(N+2); yP[1:-1] = ef[:,mode]
xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
axes[R,C].plot(xP,yP, 'blue')
axes[R,C].plot(xP,-yP, 'blue')
axes[R,C].fill_between(xP, yP,-yP,color = [1,0,1],alpha=0.2)

R = 0; C = 1; mode = 4   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_3$ = %2.2f m' % wL[mode], fontsize = 10)
yP = zeros(N+2); yP[1:-1] = ef[:,mode]
xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
axes[R,C].plot(xP,yP, 'blue')
axes[R,C].plot(xP,-yP, 'blue')
axes[R,C].fill_between(xP, yP,-yP,color = [1,0,1],alpha=0.2)

R = 1; C = 0; mode = 3   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_5$ = %2.2f m' % wL[mode], fontsize = 10)
yP = zeros(N+2); yP[1:-1] = ef[:,mode]
xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
axes[R,C].plot(xP,yP, 'blue')
axes[R,C].plot(xP,-yP, 'blue')
axes[R,C].fill_between(xP, yP,-yP,color = [1,0,1],alpha=0.2)

R = 1; C = 1; mode = 2   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_7$ = %2.2f m' % wL[mode], fontsize = 10)
yP = zeros(N+2); yP[1:-1] = ef[:,mode]
xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
axes[R,C].plot(xP,yP, 'blue')
axes[R,C].plot(xP,-yP, 'blue')
axes[R,C].fill_between(xP, yP,-yP,color = [1,0,1],alpha=0.2)

R = 2; C = 0; mode = 1   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_9$ = %2.2f m' % wL[mode], fontsize = 10)
yP = zeros(N+2); yP[1:-1] = ef[:,mode]
xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
axes[R,C].plot(xP,yP, 'blue')
axes[R,C].plot(xP,-yP, 'blue')
axes[R,C].fill_between(xP, yP,-yP,color = [1,0,1],alpha=0.2)

R = 2; C = 1; mode = 0   # x vs y 
axes[R,C].yaxis.set_tick_params(labelleft=False)
axes[R,C].set_xlim([0, L])
axes[R,C].set_ylim([-1.1, 1.1])
axes[R,C].xaxis.grid()
axes[R,C].yaxis.grid()
axes[R,C].set_title('$\lambda_{11}$ = %2.2f m' % wL[mode], fontsize = 10)
yP = zeros(N+2); yP[1:-1] = ef[:,mode]
xP = zeros(N+2); xP[-1] = L; xP[1:-1] = x;  
axes[R,C].plot(xP,yP, 'blue')
axes[R,C].plot(xP,-yP, 'blue')
axes[R,C].fill_between(xP, yP,-yP,color = [1,0,1],alpha=0.2)

fig1.savefig('a2.png')

#%% ANIMATION
nM = 6; s = 2*nM-1
mode = -1*nM + 6
yx = ef[:,mode]
Tmode = T[mode]
nT = 100
t = linspace(0,2*Tmode,nT)
yt = cos(w[mode]*t)


plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2)
fig, ax = plt.subplots(1)
#fig.subplots_adjust(top = 0.92, bottom = 0.20, left = 0.20,\
#                    right = 0.92, hspace = 0.20,wspace=0.2)

ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('y [ a.u. ]',color= 'black',fontsize = 12)
ax.set_xlabel('x  [ nm ] ',color = 'black')
ax.set_xlim([0,L])
ax.set_ylim([-1.20, 1.2])
#ax.set_xticks(np.arange(0,101,20))
#ax.set_yticks(np.arange(-20,81,20))
#txt.set_text('time = %.2f' % t[n] + ' s')  
ax.set_title('mode n = %2.0f' % s, fontsize = 12)

line1, = ax.plot([], [], 'b', lw = 2)     
ax.plot(x, yx, 'r', lw = 2)     
fig.tight_layout()

#ax.grid('visible')

def init():  
   line1.set_data([], [])
   return  line1,
   
def animate(n):
     u1 = x
     v1 = yx*yt[n]
     line1.set_data([u1], [v1])   
     
     time.sleep(0.1)
     return   line1,

anim = FuncAnimation(fig, animate, init_func = init, 
                     frames = 100, interval = 3, blit = True, repeat = False)

anim.save('agA.gif', fps = 5)


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %2.0f s' % tExe)

