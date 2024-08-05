# -*- coding: utf-8 -*-
"""
qmS02.py    Aug 2024

QUANTUM MECHANICS
Finite Difference Time Development Method
[1D] Schrodinger Equation 
   Motion of electron in an electric field   


Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH MATLAB 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmS02.pdf

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

N = 701            #  number of grid points [701]
Nt = 6000           #  Number of time steps [1 - 6000]
L = 40e-9          #  Width of X domain [40e-9] 

mEff = 1.08         # Effective mass of semiconductor
me   = 9.10938291e-31     # electron mass
hbar = 1.054571726e-34    # hbar Planck's constant
h    = 6.626e-34
e    = 1.602176565e-19    # elementary charge
m = mEff*me               # mass of wavepacket

# Scaling factors
sx = 1e-9           # position    m <---> nm
se = e              # energy      J <---. eV
st = 1e-15          # time        s <---. fs

# X grid, time steps and model parameters
x = np.linspace(0,L,N)
dx = x[2] - x[1]

c1 = 0.2
dt = c1 * 2 * m * dx**2 / hbar
tMax = Nt*dt
t =  np.arange(0,tMax,dt)          
c2 = e*dt/hbar

# Potential energy function
U0 = 2
U = zeros(N)              
U = U0*(1 - x/L)


#%% INITIAL WAVE PACKET
nx1 = round(N/4)         # pulse centre   [round(Nx/2)]
xC = x[nx1]
s = 2e-9           # pulse width    
wL = 2e-9          # wavelength

k1 = -0.5*((x-xC)/s)**2; k2 = 2*pi*(x-xC)/wL  
yR = exp(k1)#*cos(k2)
yI = exp(k1)#*sin(k2)

M = yR**2 + yI**2
A = simps(M,x)

yR = yR/sqrt(A); yI = yI/sqrt(A)

pdI = yR**2 + yI**2

R1 = yR; I1 = yI


#%%  FDTD method: Time evolution of wavefunction
def fI(I1,R1):
  I2 = zeros(N)
  for x in range(2,N-1,1):
      I2[x] = I1[x] + c1*(R1[x+1] - 2*R1[x] + R1[x-1]) - c2*U[x]*R1[x]  
  return I2 

def fR(I1,R1):
  R2 = zeros(N)
  for x in range(2,N-1,1):
      R2[x] = R1[x] - c1*(I1[x+1] - 2*I1[x] + I1[x-1]) + c2*U[x]*I1[x]  
  return R2 
 

#%%   SOLVE [2D] Schrodinger equation   
for c in range(1,Nt,1):
# Update real part of wavefunction    
  R2 = fR(I1,R1)
  R1 = R2

# Update imaginary part of wavefunction
  I2 = fI(I1,R1)

# Probability Density Function  
  probD = R2**2 + I1*I2
  I1 = I2
  
R2 = R1
I2 = (I1+I2)/2


#%% Expectation values
# Probability
w = R1 + I1*1j 
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
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.4)
fig1, ax = plt.subplots(1)
#fig.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.20,\
#                    right = 0.92, hspace = 0.20,wspace=0.2)

xP = x/sx; yP = R2/max(R2)
ax.plot(xP,yP,'b',lw = 2)
yP = I2/max(I2)
ax.plot(xP,yP,'r',lw = 1)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([-1.1, 1.1])
ax.set_ylabel('probD [ 1/m ]',color= 'black',fontsize=12)
ax.set_xlabel('x  [ nm ]',color = 'black')

ax2 = ax.twinx()
xP = x/sx; yP = U
ax2.plot(xP,yP,'k',lw = 1)
#ax.set_xlim([0, 8e-24])
ax2.set_ylim([0, 0.25])
ax2.set_ylabel('U [ eV ]',color= 'black',fontsize=11)

q = Nt*dt/st
ax2.set_title(r't = %0.f fs' %q + '   <U> = %0.3f' %Uavg
             + '  <K> = %0.3f' %Kavg + '  <E> = %0.3f  eV' %Eavg,fontsize = 11 )
fig1.tight_layout()

#%%
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,2.6)
fig2, ax = plt.subplots(1)
#fig.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.20,\
#                    right = 0.92, hspace = 0.20,wspace=0.2)

xP = x/sx; yP = probD
ax.plot(xP,yP,'b',lw = 2)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylim([0, 3e8])
ax.set_ylabel('probD [ 1/m ]',color= 'black',fontsize=12)
ax.set_xlabel('x  [ nm ]',color = 'black')

ax2 = ax.twinx()
xP = x/sx; yP = U
ax2.plot(xP,yP,'k',lw = 1)
#ax.set_xlim([0, 8e-24])
ax2.set_ylim([0, 2.1])
ax2.set_ylabel('U [ eV ]',color= 'black',fontsize=11)

q = Nt*dt/st
ax2.set_title(r't = %0.f fs' %q + '   <U> = %0.4f' %Uavg
             + '  <K> = %0.4f' %Kavg + '  <E> = %0.4f  eV' %Eavg,fontsize = 11 )
fig2.tight_layout()


#%% SAVE FIGURES
fig1.savefig('a1.png')
fig2.savefig('a2.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)







