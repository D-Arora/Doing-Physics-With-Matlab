# -*- coding: utf-8 -*-
"""
qm0020.py    April 2024

QUANTUM MECHANICS
Finite Difference Time Development Method: Animation
     [1D] Schrodinger Equation:
     Free particle: Scattering and Tunnelling

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH PYTHON 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm020.htm

"""

#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag
from scipy.integrate import odeint, quad, dblquad, simps

tStart = time.time()


#%% INPUT VARIABLES
Nx = 701             #  ODD number - number of grid points [701]
Nt = 1200            #  Number of time steps [30000]
L = 40e-9             #  Half-width of X domain [5e-9 nm] 
f = 100              # number of frames

U0 = 60              # Potential energy  [eV]

C1 = 1/5             # Schrodinger equation constant

# INITIAL WAVE PACKET
s = 2e-10               # pulse width  [m]
wL = 1.5e-10            # nominal wavelength  [  1.5e-10 m]
xC = -4e-9              # pulse centre [m]



#%% CONSTANTS AND VARIABLES
me   = 9.10938291e-31     # electron mass
hbar = 1.054571726e-34  # hbar Planck's constant
h    = 6.626e-34
e    = 1.602176565e-19     # elementary charge
k1 = -hbar**2/(2*me)

#% Initialize arrays  
# Wavefunction y yR yI / prob density pd  
psiR = zeros([Nt,Nx])      # real part wavefunction
psiI = zeros([Nt,Nx])      # imaginary part wvefunction
psi = psiR + 1j*psiI       # complex wavefunction
pd = zeros([Nt,Nx])        # probability density
Uavg = zeros([Nt])         # expectation value: potential energy  [J]
Kavg = zeros([Nt])         # expectation value: kinetic energy    [J]
Eavg = zeros([Nt])         # expectation value: total energy      [J]
Prob = zeros([Nt])         # Probability


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


#%% SETUP CALCULATIONS
x = np.linspace(-L,L,Nx); dx = x[2] - x[1]
dt = C1 * 2 * me * dx**2 / hbar
tMax = Nt*dt
t =  np.arange(0,tMax,dt)           # 0:dt:Nt*dt;
C2 = e*dt/hbar;
C3 = -hbar**2 / (2 * me * dx**2 * e)

# POTENTIAL: zero, step, barrier at x = 0   [eV]
U = zeros(Nx)              
U[x>0] = U0
#U[x>0.2e-9] = 0       # 0.2e-9 m

# INITIAL WAVEPACKET and NORMALIZE WAVEPACKET
k1 = -0.5*((x-xC)/s)**2; k2 = 2*pi*(x-xC)/wL  
yR = exp(k1)*cos(k2)
yI = exp(k1)*sin(k2)

fn = yR**2 + yI**2
A = simps(fn,x)

yR = yR/sqrt(A); yI = yI/sqrt(A)    # Initial normalized wavepacket

pdI = yR**2 + yI**2       # Initial probability density 


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

#Expectation values and probabilities
for nt in range(0,Nt):
    w = psiR[nt,:] + psiI[nt,:]*1j 
    
    fn = np.real(np.conj(w)*w)
    Prob[nt] = simps(fn,x)
    
    y2 = secondDer(Nx,dx)@w
    fn = np.real(np.conj(w)*y2)
    Kavg[nt] = -hbar**2*simps(fn,x)/(2*me)
 
    fn = np.real(np.conj(w)*U*w)
    Uavg[nt] = simps(fn,x)
    
Eavg = Kavg/e + Uavg

# Percentage probailities
fn = np.real(np.conj(w)*w)
Prob0 = simps(fn,x)
w1 = w[x<0]; w2 = w[x>0]
fn = np.real(np.conj(w1)*w1)
Prob1 = simps(fn,x[x<0])*100
fn = np.real(np.conj(w2)*w2)
Prob2 = simps(fn,x[x>0])*100
prob = Prob1 + Prob2

#%% CONSOLE OUTPUT
print('  ')
print('Prob( x < 0) = %2.1f' %Prob1)
print('Prob( x > 0) = %2.1f' %Prob2)
print('Prob         = %2.1f' %prob)
print('  ')


#%%  Time evolution plots
plt.rcParams["figure.figsize"] = (4,4)
fig1, axes = plt.subplots(nrows=1, ncols=1)
fig1.subplots_adjust(top = 0.94, bottom = 0.15, left = 0.180,\
                    right = 0.92, hspace = 0.36,wspace=0.40)

tP = t*1e15      # [fs]
    
# t vs energy
axes.set_xlabel('t  [ fs ]',color= 'black',fontsize = 12)
axes.set_ylabel('< energy >   [ eV ] ',color = 'black',fontsize = 12)
#axes[R].set_xlim([0, 0.7])
#axes[R].set_ylim([0, 5])
#axes[R,C].set_xticks(np.arange(0,101,20))
#axes[R,C].set_yticks(np.arange(-20,81,20))
axes.xaxis.grid()
axes.yaxis.grid()
yP = Kavg/e
axes.plot(tP,yP,'b',lw = 2, label = 'K')
yP = Uavg
axes.plot(tP,yP,'r',lw = 2, label = 'U')
yP = Kavg/e + Uavg
axes.plot(tP,yP,'k',lw = 2, label = 'E')
axes.legend()

fig1.savefig('a1.png')



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
ax2.set_ylim([0, 100])
#ax2.set_ylim([-100, 0])
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

anim.save('ag.gif', fps = 6)  


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)