# -*- coding: utf-8 -*-
"""
qmS01.py    July 024

QUANTUM MECHANICS
Finite Difference Time Development Method
[1D] Schrodinger Equation 
   Conduction Band of a semiconductor   


Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH MATLAB 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmS01.pdf

"""


#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag
from scipy.integrate import odeint, quad, dblquad, simps

tStart = time.time()

#%% FUNCTIONS
def firstDer(N,dx):
    v  = ones(N-1)
    M1 = diag(-v,-1)
    M2 = diag(v,1)
    M = M1+M2
    M[0,0] = -2; M[0,1] = 2; M[N-1,N-2] = -2; M[N-1,N-1] = 2
    MF = M/(2*dx) 
    return MF

def secondDer(N,dx):
     v = -2*ones(N)
     M1 = np.diag(v)
     v = np.ones(N-1)
     M2 = np.diag(v,1)
     M3 = np.diag(v,-1)
     M = M1+M2+M3
     M[0,0] = 1; M[0,1] = -2; M[0,2] = 1
     M[N-1,N-3] = 1; M[N-1,N-2] = -2; M[N-1,N-1]=1
     MS = M/(dx**2) 
     return MS


#%% SETUP

N = 999            #  number of grid points [701]
Nt = 20000      #  Number of time steps [40]
L = 40e-9           #  Width of X domain [40e-9] 

mEff = 1.08         # Effective mass of semiconductor
me   = 9.10938291e-31     # electron mass
hbar = 1.054571726e-34  # hbar Planck's constant
h    = 6.626e-34
e    = 1.602176565e-19     # elementary charge
m = mEff*me

# Scaling factors
sx = 1e-9           # position    m <---> nm
se = 1e-19          # energy      J <---. eV
st = 1e-15          # time        s <---. fs

x = np.linspace(0,L,N); dx = x[2] - x[1]
k1 = -hbar**2/(2*me)
C1 = 0.1
dt = C1 * 2 * m * dx**2 / hbar
tMax = Nt*dt
t =  np.arange(0,tMax,dt)           # 0:dt:Nt*dt
C2 = e*dt/hbar;

pd = zeros(N)            # probability density


U = zeros(N)          # potential energy
nU = round(N/2)
U[0:nU] = 0.1
U[nU:N] = 0.2

# Prob = zeros([Nt])       # probability
# Uavg = zeros([Nt])       # expectation value: potential energy
# Kavg = zeros([Nt])       # expectation value: kinetic energy
# Eavg = zeros([Nt])       # total energy
# xavg = zeros([Nt])       # expectation value: position
# pavg = zeros([Nt])       # momentum
# vavg = zeros([Nt])       # velocity
# deltaX = zeros([Nt])
# delta = zeros([Nt])

# INITIAL WAVE PACKET
nx1 = round(N/4)         # pulse centre   [round(Nx/2)]
s = 2e-9           # pulse width    [L/20]
wL = 2e-9               # wavelength

k1 = -0.5*((x-x[nx1])/s)**2; k2 = 2*pi*(x-x[nx1])/wL  
yR = exp(k1)*cos(k2)
yI = exp(k1)*sin(k2)

M = yR**2 + yI**2
z = linspace(0,L,N)
A = simps(M,z)

yR = yR/sqrt(A); yI = yI/sqrt(A)

pdI = yR**2 + yI**2


#%%
# Solve Schrodinger Equation: FDTD Method

for nt in range(Nt):
   for nx in range(1,N-2): 
       yR[nx] = yR[nx] - C1*(yI[nx+1]-2*yI[nx]+yI[nx-1]) + C2*U[nx]*yI[nx]
      # psiR[nt,:] = yR
      
   for nx in range(1,N-2):
      yI[nx] = yI[nx] + C1*(yR[nx+1]-2*yR[nx]+yR[nx-1]) - C2*U[nx]*yR[nx]
    
pd = yR**2 + yI**2


#%% Expectation values
# Probability
w = yR + yI*1j 
fn = np.real(np.conj(w)*w)
Prob = simps(fn,x)
# Potential energy
fn = np. conj(w)*U*w
Uavg = np.real(simps(fn,x))
# Kinetic energy
y2 = secondDer(N,dx)@w
fn = np.conj(w)*y2
Kavg = np.real(-hbar**2*simps(fn,x)/(2*m*se))
# Total energy
Eavg = Kavg + Uavg

# CONSOLE OUTPUT
print('Expectation values')
print(r'   prob = %2.3f   eV ' % Prob)
q = Nt*dt/st; print(r'   t = %2.3f  fs ' % q)
print(r'   <U> = %2.3f   eV ' % Uavg)
print(r'   <K> = %2.3f   eV ' % Kavg)
print(r'   <E> = %2.3f   eV ' % Eavg)
q = (h/wL)**2/(2*m*se); print(r'   KEpulse = %2.3f   eV ' % q)


#%% GRAPHICS
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(1)
#fig.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.20,\
#                    right = 0.92, hspace = 0.20,wspace=0.2)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('probD [ 1/m ]',color= 'black')
ax.set_xlabel('x  [ nm ]',color = 'black')
ax.plot(x/sx,pd,'b',lw = 2)
#ax.plot(x/sx,yI,'r',lw = 2)

fig1.tight_layout()
fig1.savefig('a1.png')
#ax.set_xlim([0, 8e-24])
#ax.set_ylim([-20, 80])
#ax.set_xticks(np.arange(0,101,20))
#ax.set_yticks(np.arange(-20,81,20))
#ax.set_title('Trajectory', fontsize = 12)
#ax.text(1, 70, 'y1$_{max}$ = %2.2f m' % max(y1), fontsize = 12, color = 'blue')
#ax.text(1, -10, 'y2$_{max}$ = %2.2f m' % max(y2), fontsize = 12, color = 'red')
# fig.tight_layout()
# xP = h*K; yP = psd
# ax.plot(xP,yP,'b',lw = 2)
# plt.rcParams["figure.figsize"] = (7,7)
# fig1, axes = plt.subplots(nrows=2, ncols=2)
# fig1.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.120,\
#                     right = 0.95, hspace = 0.36,wspace=0.40)

# tP = t*1e15
    
# R = 0; C = 0   # t vs K
# axes[R,C].set_xlabel('t  [ fs ]',color= 'black',fontsize = 12)
# axes[R,C].set_ylabel('< K >  [ eV ]  ',color = 'black',fontsize = 12)
# axes[R,C].set_xlim([0, 0.7])
# axes[R,C].set_ylim([0, 70])
# #axes[R,C].set_xticks(np.arange(0,101,20))
# #axes[R,C].set_yticks(np.arange(-20,81,20))
# axes[R,C].xaxis.grid()
# axes[R,C].yaxis.grid()
# yP = Kavg/e
# axes[R,C].plot(tP,yP,'b',lw = 2)

# R = 0; C = 1   # t vs x,y
# axes[R,C].set_xlabel('t  [ fs ]',color= 'black',fontsize = 12)
# axes[R,C].set_ylabel('< x >   [ nm ] ',color = 'black',fontsize = 12)
# axes[R,C].set_xlim([0, 0.7])
# axes[R,C].set_ylim([0, 5])
# #axes[R,C].set_xticks(np.arange(0,101,20))
# #axes[R,C].set_yticks(np.arange(-20,81,20))
# axes[R,C].xaxis.grid()
# axes[R,C].yaxis.grid()
# yP = xavg*1e9
# axes[R,C].plot(tP,yP,'b',lw = 2)

# R = 1; C = 0  # t vs <p>
# axes[R,C].set_ylabel('< p >  [ x10$^{23}$  N.s ]',color= 'black',fontsize = 12)
# axes[R,C].set_xlabel('t  [ fs ]',color = 'black',fontsize = 12)
# axes[R,C].set_xlim([0, 0.7])
# axes[R,C].set_ylim([0, 1.0])
# # axes[R,C].set_xticks(np.arange(0,11,2))
# axes[R,C].xaxis.grid()
# axes[R,C].yaxis.grid()
# yP = sqrt(2*me*Kavg)*1e23
# axes[R,C].plot(tP,yP, 'b',lw = 2)

# R = 1; C = 1  # t vs delta
# axes[R,C].set_ylabel('delta/hbar',color= 'black',fontsize = 12)
# axes[R,C].set_xlabel('t  [ fs ]',color = 'black',fontsize = 12)
# axes[R,C].set_xlim([0, 0.7])
# # axes[R,C].set_xticks(np.arange(0,11,2))
# axes[R,C].xaxis.grid()
# axes[R,C].yaxis.grid()
# yP = delta/hbar
# axes[R,C].plot(tP,yP, 'b',lw = 2)
# axes[R,C].plot([0,np.amax(tP)],[0.5,0.5], 'r',lw = 1)

# fig1.savefig('a1.png')


# #%%   Fourier transform PSI at t = 0   K = 1/wL
# KMax = 2/wL; KMin= 0; nK = 2001
# K = linspace(KMin,KMax,nK);
# hP = psiR[-10,:]   # pd[-1000,:] 
# HP = zeros(nK);HPR = zeros(nK); HPI = zeros(nK)
# for c in range(nK-1):
#      g = hP* exp(1j*2*pi*K[c]*x)
#      gR = np.real(g); gI = np.imag(g)
#      HPR[c] = simps(gR,x); HPI[c] = simps(gI,x)
#     # HP[c] = simps(g,x)

# HP = HPR + 1j*HPI
# #psd = HPI*HPR
# #psd = psd/max(psd)
# psd = np.conj(HP)*HP
# psd = psd/max(psd)



# # fig.savefig('a2.png')


# #%% GRAPHICS: INITIALIZE ANIMATION PLOTS 

# xP = x*1e9

# pdMax = np.amax(pd)*1e-9

# yRmax = 1.1*np.amax(psiR); yRmin = 1.1*np.amin(psiR)
# yImax = 1.1*np.amax(psiI); yImin = 1.1*np.amin(psiI)

# plt.rcParams['font.family'] = ['Tahoma']
# plt.rcParams['font.size'] = 12

# plt.rcParams["figure.figsize"] = (5,8)
# fig, ax = plt.subplots(nrows=3, ncols=1)
# fig.subplots_adjust(top = 0.97, bottom = 0.08, left = 0.22,\
#                     right = 0.93, hspace = 0.36,wspace=0.40)

# R = 0  # x vs y 
# ax[R].set_xlabel('$x$  [nm]',color= 'black',fontsize = 12)
# ax[R].set_ylabel('$\psi_R$',color = 'black',fontsize = 14)
# ax[R].set_xlim([0, L*1e9])
# ax[R].set_ylim([yRmin,yRmax])
# ax[R].xaxis.grid()
# ax[R].yaxis.grid()
# line0, = ax[R].plot([], [], 'b', lw = 1)
# ax[R].plot(xP, psiR[0,:],'r', lw = 1) 
# ax[R].grid('visible') 

# R = 1  # t vs deltaX 
# ax[R].set_xlabel('$x $ [nm]',color= 'black',fontsize = 12)
# ax[R].set_ylabel('$\psi_I$',color = 'black',fontsize = 14)
# ax[R].set_xlim([0, L*1e9])
# ax[R].set_ylim([yImin,yImax])
# #ax[R].set_xticks(np.arange(0,101,20))
# #ax[R].set_yticks(np.arange(-20,81,20))
# ax[R].xaxis.grid()
# ax[R].yaxis.grid()
# line1, = ax[R].plot([], [], 'b', lw = 1)
# ax[R].plot(xP, psiI[0,:],'r', lw = 1) 
# ax[R].grid('visible') 

# R = 2  # x vs y 
# ax[R].set_xlabel('$x$  [nm]',color= 'black',fontsize = 12)
# ax[R].set_ylabel('$|\psi|^2$  [$nm^{-1}$]',color = 'black',fontsize = 12)
# ax[R].set_xlim([0, L*1e9])
# ax[R].set_ylim([0,pdMax])
# #ax[R].set_xticks(np.arange(0,101,20))
# ax[R].set_yticks(np.arange(0,3,1))
# ax[R].xaxis.grid()
# ax[R].yaxis.grid()
# line, = ax[R].plot([], [], 'b', lw = 2)
# ax[R].plot(xP, pdI*1e-9,'r', lw = 1) 
# ax[R].grid('visible') 
# time_text = ax[R].text(8e-9,0.8*pdMax, '')


# #%% FUNCTIONS
# def init():  
#     line.set_data([], [])
    
#     # time_text.set_text('')
    
#     return line,
   
# def animate(n):
#       m = round(n*Nt/f)
#       u = xP
      
#       v = pd[m,:]*1e-9
#       line.set_data([u], [v]) 
    
#       w = psiR[m,:]
#       line0.set_data([u], [w]) 
      
#       z = psiI[m,:]
#       line1.set_data([u], [z]) 
      
    
#       T = t[m]*1e15
#       time_text.set_text('   time = %.3f' % T + ' fs')  
    
#       time.sleep(0.1)
#       return   line, line0, line1,  time_text,

    
# anim = FuncAnimation(fig, animate, init_func = init, 
#                       frames = f, interval = 2, blit = True, repeat = False)

# #anim.save('ag.gif', fps = 12)  


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)







