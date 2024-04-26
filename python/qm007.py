# -*- coding: utf-8 -*-
"""
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

Nx = 701             #  ODD number - number of grid points [701]
Nt = 8200            #  Number of time steps [30000]
L = 5e-9             #  Width of X domain [5e-9] 
f = 100              # number of frames

me   = 9.10938291e-31     # electron mass
hbar = 1.054571726e-34  # hbar Planck's constant
h    = 6.626e-34
e    = 1.602176565e-19     # elementary charge

x = np.linspace(-L,L,Nx); dx = x[2] - x[1]
k1 = -hbar**2/(2*me)
C1 = 1/5;
dt = C1 * 2 * me * dx**2 / hbar
tMax = Nt*dt
t =  np.arange(0,tMax,dt)           # 0:dt:Nt*dt;
C2 = e*dt/hbar;
C3 = -hbar**2 / (2 * me * dx**2 * e)

#% Initialize arrays  
# Wavefunction y yR yI / prob density pd  
psiR = zeros([Nt,Nx])      # real part wavefunction
psiI = zeros([Nt,Nx])      # imaginary part wvefunction
psi = psiR + 1j*psiI       # complex wavefunction
pd = zeros([Nt,Nx])        # probability density
phase = zeros([Nt,Nx])     # phase
#z = zeros([Nt,Nx])

U = zeros(Nx)          # potential energy  [eV]
Uavg = zeros([Nt])       # expectation value: potential energy
Kavg = zeros([Nt])       # expectation value: kinetic energy
Eavg = zeros([Nt])       # total energy
xavg = zeros([Nt])       # expectation value: position
pavg = zeros([Nt])       # momentum
vavg = zeros([Nt])       # velocity
deltaX = zeros([Nt])
delta = zeros([Nt])
Prob = zeros([Nt])

# potential energy  [eV]
U[x>0] = 70
a = 6.25*1e18
U = 1*a*x**2

# INITIAL WAVE PACKET
nx1 = round(Nx/10)        # pulse centre   [round(Nx/2)]
s = 2e-10   #L/25 #L/15              # pulse width    [L/20]
wL = 1.5e-10          # wavelength
xC = 0   #x[nx1]
k1 = -0.5*((x-xC)/s)**2; k2 = 2*pi*(x-xC)/wL  
yR = exp(k1)*cos(k2)
yI = exp(k1)*sin(k2)

M = yR**2 + yI**2
z = linspace(0,L,Nx)
A = simps(M,z)

yR = yR/sqrt(A); yI = yI/sqrt(A)

pdI = yR**2 + yI**2


#%% Solve Schrodinger Equation: FDTD Method

for nt in range(Nt):
   for nx in range(1,Nx-2): 
       yR[nx] = yR[nx] - C1*(yI[nx+1]-2*yI[nx]+yI[nx-1]) + C2*U[nx]*yI[nx]
       psiR[nt,:] = yR
      
   for nx in range(1,Nx-2):
      yI[nx] = yI[nx] + C1*(yR[nx+1]-2*yR[nx]+yR[nx-1]) - C2*U[nx]*yR[nx]
      psiI[nt,:] = yI
      pd[nt,:] = yR**2 + yI**2
        
   
for nt in range(1,Nt-1):
   psiI[nt,:] = 0.5*(psiI[nt,:] + psiI[nt-1,:])
   pd[nt,:] = psiR[nt,:]**2 + psiI[nt,:]**2

#Expectation values and Uncertainty Principle
for nt in range(0,Nt):
    w = psiR[nt,:] + psiI[nt,:]*1j 
    
    fn = np.real(np.conj(w)*w)
    Prob[nt] = simps(fn,x)
    
    fn = np.real(np.conj(w)*x*w)
    xavg[nt] = simps(fn,x)
    
    fn = np.real(np.conj(w)*x*x*w)
    X2 = simps(fn,x)
    
    deltaX[nt] = sqrt(X2 - xavg[nt]**2)
    
    y1 = firstDer(Nx,dx)@w
    fn = np.real(np.conj(w)*y1)
    pavg = -1j*hbar*simps(fn,x)
   
    y2 = secondDer(Nx,dx)@w
    fn = np.real(np.conj(w)*y2)
    p2avg = -hbar**2*simps(fn,x)
    deltap = sqrt(p2avg - np.imag(pavg)**2)

    delta[nt] = deltaX[nt]*deltap

    Kavg[nt] = -hbar**2*simps(fn,x)/(2*me)
 
    fn = np.conj(w)*U*w
    Uavg[nt] = simps(fn,x)
    
Eavg = Kavg/e + Uavg

#Classical calculations
vavg = (xavg[-10]-xavg[10]) / (t[-10]-t[10])
Kc = 0.5*me*vavg**2/e
pc = me*vavg
vd = h/(wL*me)

# Expection values
Ke = np.amax(Kavg)/e
pe = sqrt(2*me*Ke*e)
ve = pe/me

# Group and phase velocities
vGroup = vavg
p0 = me*vGroup
f0 = Ke*e/h
vPhase = wL*f0
vR = vGroup/vPhase

#CONSOLE OUTPUT
print('Classical values')
print(r'   v = %2.2e   m/s ' % vavg)
print(r'   p = %2.2e   N.s ' % pc)
print(r'   K = %2.2f   eV ' % Kc)

print('Expectation values')
print(r'   v = %2.2e   m/s ' % ve)
print(r'   p = %2.2e   N.s ' % pe)
print(r'   K = %2.2f   eV ' % Ke)

print('Group and Phase velocities')
print(r'   wL_0 = %2.2e   m ' % wL)
print(r'   p_0  = %2.2e   N.s ' % p0)
print(r'   f_0  = %2.2e   Hz ' % f0)
print(r'   vGroup  = %2.2e   m/s ' % vGroup)
print(r'   vPhase  = %2.2e   m/s ' % vPhase)
print(r'   vGroup / vPhase  = %2.2f  ' % vR)


#%%  Time evolution plots
plt.rcParams["figure.figsize"] = (5,5)
fig1, axes = plt.subplots(nrows=2, ncols=1)
fig1.subplots_adjust(top = 0.94, bottom = 0.15, left = 0.180,\
                    right = 0.92, hspace = 0.36,wspace=0.40)

tP = t*1e15
    
R = 0;    # t vs xavg
axes[R].set_xlabel('t  [ fs ]',color= 'black',fontsize = 12)
axes[R].set_ylabel('< x >  [ nm ]  ',color = 'black',fontsize = 12)
#axes[R].set_xlim([0, 0.7])
#axes[R].set_ylim([0, 70])
#axes[R,C].set_xticks(np.arange(0,101,20))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
yP = xavg*1e9
axes[R].plot(tP,yP,'b',lw = 2)

R = 1   # t vs energy
axes[R].set_xlabel('t  [ fs ]',color= 'black',fontsize = 12)
axes[R].set_ylabel('< energy >   [ eV ] ',color = 'black',fontsize = 12)
#axes[R].set_xlim([0, 0.7])
#axes[R].set_ylim([0, 5])
#axes[R,C].set_xticks(np.arange(0,101,20))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes[R].xaxis.grid()
axes[R].yaxis.grid()
yP = Kavg/e
axes[R].plot(tP,yP,'b',lw = 2, label = 'K')
yP = Uavg
axes[R].plot(tP,yP,'r',lw = 2, label = 'U')
yP = Kavg/e + Uavg
axes[R].plot(tP,yP,'k',lw = 2, label = 'E')
axes[R].legend()

fig1.savefig('a1.png')


#%%
# #%%   Fourier transform PSI at t = 0   K = 1/wL
# KMax = 2/wL; KMin= 0; nK = 2001
# K = linspace(KMin,KMax,nK);
# hP = psiR[-10,:]   # pd[-1000,:] 
# HP = zeros(nK);HPR = zeros(nK); HPI = zeros(nK)
# for c in range(nK-1):
#       g = hP* exp(1j*2*pi*K[c]*x)
#       gR = np.real(g); gI = np.imag(g)
#       HPR[c] = simps(gR,x); HPI[c] = simps(gI,x)
#     # HP[c] = simps(g,x)

# HP = HPR + 1j*HPI
# #psd = HPI*HPR
# #psd = psd/max(psd)
# psd = np.conj(HP)*HP
# psd = psd/max(psd)

# plt.rcParams["figure.figsize"] = (5,3)
# fig, ax = plt.subplots(1)
# #fig.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.20,\
# #                    right = 0.92, hspace = 0.20,wspace=0.2)
# ax.xaxis.grid()
# ax.yaxis.grid()
# ax.set_ylabel('psd  [ a. u. ]',color= 'black')
# ax.set_xlabel('p  [ N.s ]',color = 'black')
# ax.set_xlim([0, 8e-24])
# #ax.set_ylim([-20, 80])
# #ax.set_xticks(np.arange(0,101,20))
# #ax.set_yticks(np.arange(-20,81,20))
# #ax.set_title('Trajectory', fontsize = 12)
# #ax.text(1, 70, 'y1$_{max}$ = %2.2f m' % max(y1), fontsize = 12, color = 'blue')
# #ax.text(1, -10, 'y2$_{max}$ = %2.2f m' % max(y2), fontsize = 12, color = 'red')
# fig.tight_layout()
# xP = h*K; yP = psd
# ax.plot(xP,yP,'b',lw = 2)
# fig.savefig('a2.png')


#%% GRAPHICS: INITIALIZE ANIMATION PLOTS 

xP = x*1e9

pdMax = np.amax(pd)*1e-9

yRmax = 1.1*np.amax(psiR); yRmin = 1.1*np.amin(psiR)
yImax = 1.1*np.amax(psiI); yImin = 1.1*np.amin(psiI)

plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

plt.rcParams["figure.figsize"] = (5,7)
fig, ax = plt.subplots(nrows=3, ncols=1)
fig.subplots_adjust(top = 0.97, bottom = 0.08, left = 0.22,\
                    right = 0.84, hspace = 0.36,wspace=0.40)

R = 0  # x vs psiR
yP = psiR[0,:] / np.amax(psiR)
ax[R].set_xlabel('$x$  [nm]',color= 'black',fontsize = 12)
ax[R].set_ylabel('$\psi_R$',color = 'black',fontsize = 14)
ax[R].set_xlim([-L*1e9, L*1e9])
ax[R].set_ylim([-1.1,1.1])
ax[R].xaxis.grid()
ax[R].yaxis.grid()
line0, = ax[R].plot([], [], 'b', lw = 1)
ax[R].plot(xP, yP,'r', lw = 1) 
ax[R].grid('visible') 

R = 1  # t vs psiI
yP = psiR[0,:] / np.amax(psiR)
ax[R].set_xlabel('$x $ [nm]',color= 'black',fontsize = 12)
ax[R].set_ylabel('$\psi_I$',color = 'black',fontsize = 14)
ax[R].set_xlim([-L*1e9, L*1e9])
ax[R].set_ylim([-1.1,1.1])
#ax[R].set_xticks(np.arange(0,101,20))
#ax[R].set_yticks(np.arange(-20,81,20))
ax[R].xaxis.grid()
ax[R].yaxis.grid()
line1, = ax[R].plot([], [], 'b', lw = 1)
ax[R].plot(xP, yP,'r', lw = 1) 
ax[R].grid('visible') 

R = 2  # x vs y 
ax[R].set_xlabel('$x$  [nm]',color= 'black',fontsize = 12)
ax[R].set_ylabel('$|\psi|^2$  [$nm^{-1}$]',color = 'black',fontsize = 12)
ax[R].set_xlim([-L*1e9, L*1e9])
ax[R].set_ylim([0,pdMax])
#ax[R].set_xticks(np.arange(0,101,20))
#ax[R].set_yticks(np.arange(0,3,1))
ax[R].xaxis.grid()
ax[R].yaxis.grid()
line, = ax[R].plot([], [], 'b', lw = 2)
ax[R].plot(xP, pdI*1e-9,'r', lw = 1) 
ax[R].grid('visible') 
time_text = ax[R].text(0,0.8*pdMax, '')

ax2 = ax[R].twinx()
ax2.set_ylabel('U   [eV]',color= 'black')

#ax2.set_xlim([0, 11])
ax2.set_ylim([0, 150])
#ax2.set_xticks(np.arange(0,11,2))
#ax2.set_yticks(np.arange(-80,61,20))
#ax2.tick_params(axis='y', labelcolor = 'red')
#ax2.yaxis.grid(color = 'r')
ax2.plot(x*1e9,U,'k',lw = 1)


#%% FUNCTIONS
def init():  
    line.set_data([], [])
    
    # time_text.set_text('')
    
    return line,
   
def animate(n):
      m = round(n*Nt/f)
      u = xP
      
      v = pd[m,:]*1e-9 
      line.set_data([u], [v]) 
    
      w = psiR[m,:] / np.amax(psiR)
      line0.set_data([u], [w]) 
      
      z = psiI[m,:] / np.amax(psiI)
      line1.set_data([u], [z]) 
      
    
      T = t[m]*1e15
      time_text.set_text('   time = %.3f' % T + ' fs')  
    
      time.sleep(0.1)
      
      return   line, line0, line1,  time_text,

#%%
anim = FuncAnimation(fig, animate, init_func = init, 
                      frames = f, interval = 2, blit = True, repeat = False)

#%%
anim.save('ag.gif', fps = 12)  



#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)