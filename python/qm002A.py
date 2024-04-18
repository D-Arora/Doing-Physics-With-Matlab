# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 18:11:28 2024

@author: Owner


# qm002.py    April 2024


"""

#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt
from scipy.integrate import odeint, quad, dblquad, simps

tStart = time.time()

Nx = 199                 #  ODD number - number of grid points [701]
Nt = 2000              #  Number of time steps [30000]
L = 100e-9                 #  Width of X domain [5e-9] 
f = 200                 # number of frames

me = 9.10938291e-31     # electron mass
hbar = 1.054571726e-34  # hbar Planck's constant
e = 1.602176565e-19     # elementary charge

x = np.linspace(0,L,Nx); dx = x[2] - x[1]
k1 = -hbar**2/(2*me)
C1 = 1/10;
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
z = zeros([Nt,Nx])

U = zeros(Nx)          # potential energy
Prob = zeros([1,Nt])       # probability
Uavg = zeros([1,Nt])       # expectation value: potential energy
Kavg = zeros([1,Nt])       # expectation value: kinetic energy
Eavg = zeros([1,Nt])       # total energy
xavg = zeros([1,Nt])       # expectation value: position
pavg = zeros([1,Nt])       # momentum
vavg = zeros([1,Nt])       # velocity


# INITIAL WAVE PACKET

nx1 = round(Nx/2)        # pulse centre   [round(Nx/2)]
s = 5e-9 #L/15              # pulse width    [L/20]
wL = 1.6e-10          # wavelength
  
yR = exp(-0.5*((x-x[nx1])/s)**2)
#yI = exp(-0.5.*((x-x(nx1))./s).^2)
yI = zeros(Nx)

M = yR**2 + yI**2
z = linspace(0,L,Nx)
A = simps(M,z)

yR = yR/sqrt(A); yI = yI/sqrt(A)

pdI = yR**2 + yI**2



# Solve Schrodinger Equation: FDTD Method

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

#plt.plot(x,pd[15000-1,:]) 

xP = x*1e9

pdMax = np.amax(pd)*1e-9

yRmax = 1.1*np.amax(psiR); yRmin = 1.1*np.amin(psiR)
yImax = 1.1*np.amax(psiI); yImin = 1.1*np.amin(psiI)

#%%  # Initializing a figure in which the graph will be plotted 
# font1 = {'family':'Tahoma','color':'black','size':12}
# plt.rcParams['font.family'] = ['Tahoma']
# plt.rcParams['font.size'] = 12

plt.rcParams["figure.figsize"] = (5,8)
fig, ax = plt.subplots(nrows=3, ncols=1)
fig.subplots_adjust(top = 0.94, bottom = 0.12, left = 0.180,\
                    right = 0.95, hspace = 0.36,wspace=0.40)


R = 0  # x vs y 
ax[R].set_xlabel('$x  [nm]$',color= 'black',fontsize = 12)
ax[R].set_ylabel('$\psi_R$',color = 'black',fontsize = 14)
ax[R].set_xlim([0, L*1e9])
ax[R].set_ylim([yRmin,yRmax])
ax[R].xaxis.grid()
ax[R].yaxis.grid()
line0, = ax[R].plot([], [], 'b', lw = 3)
ax[R].plot(xP, psiR[0,:],'r', lw = 1) 
ax[R].grid('visible') 


R = 1  # x vs y 
ax[R].set_xlabel('$x  [nm]$',color= 'black',fontsize = 12)
ax[R].set_ylabel('$\psi_I$',color = 'black',fontsize = 14)
ax[R].set_xlim([0, L*1e9])
ax[R].set_ylim([yImin,yImax])
#ax[R].set_xticks(np.arange(0,101,20))
#ax[R].set_yticks(np.arange(-20,81,20))
ax[R].xaxis.grid()
ax[R].yaxis.grid()
line1, = ax[R].plot([], [], 'b', lw = 3)
ax[R].plot(xP, psiI[0,:],'r', lw = 1) 
ax[R].grid('visible') 


R = 2  # x vs y 
ax[R].set_xlabel('$x  [nm]$',color= 'black',fontsize = 12)
ax[R].set_ylabel('$|\psi|^2$  [$nm^{-1}$]',color = 'black',fontsize = 12)
ax[R].set_xlim([0, L*1e9])
ax[R].set_ylim([0,pdMax])
#ax[R].set_xticks(np.arange(0,101,20))
ax[R].set_yticks(np.arange(0,0.13,0.020))
ax[R].xaxis.grid()
ax[R].yaxis.grid()
line, = ax[R].plot([], [], 'b', lw = 3)
ax[R].plot(xP, pdI*1e-9,'r', lw = 1) 
ax[R].grid('visible') 
time_text = ax[R].text(8e-9,0.8*pdMax, '')


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
    
      w = psiR[m,:]
      line0.set_data([u], [w]) 
      
      z = psiI[m,:]
      line1.set_data([u], [z]) 
      
    
      T = t[n]*1e15
      time_text.set_text('   time = %.1f' % T + ' fs')  
    
      time.sleep(0.1)
      return   line, line0, line1,  time_text,
 
anim = FuncAnimation(fig, animate, init_func = init, 
                      frames = f, interval = 2, blit = True, repeat = False)

tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

# anim.save('ag.gif', fps = 2)  





