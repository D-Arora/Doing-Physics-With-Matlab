# -*- coding: utf-8 -*-
"""
qm031.py            May 2024

Ian Cooper         matlabvisualphysics@gmail.com

Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

QUANTUM MECHANICS
   Scattering from a finite setp potential
   Analtyical solutions
   
   https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm030.pdf

"""

# LIBRARIES  ==============================================================
import numpy
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, real, imag
import cmath
import matplotlib.pyplot as plt


#%%  CONSTANTS and VARIABLES
e = 1.6021766208e-19           # Elementary charge [C]
h = 6.626070040e-34            # Planck constant [ J.s]
hbar = 1.054571800e-34         # Reduced Planck's constant
me = 9.10938356e-31            # Electron mass [kg]

xMin = -0.6e-9; xMax = 0.6e-9     # X domain  [m]
E0 = 80                   # beam particle energy  [eV]
U0 = 100                  # height of potential step [eV]
N = 8299                           # Grid points

            
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

probD1 = psi1R**2 + psi1I**2
probD2 = psi2R**2 + psi2I**2


# Velocities, fluxes, reflection coeff, transmission coeff
v1 = hbar*k1/me; v2 = hbar*k2/me
Jinc   = v1*A**2
Jref   = -v1*B**2
Jtrans = v2*C**2

Tc = real(4*k1*k2/(k1+k2)**2)
Rc = 1 - Tc


#%% CONSOLE OUTPUT
print(' ')
print(r'   E  = %2.1f  eV ' % E0)
print(r'   U0 = %2.1f  eV ' % U0)
print(r'   K1  = %2.0f  eV ' % E0)
s = E0 - U0; print(r'   K2  = %2.0f  eV ' % s)
s = E0*e/hbar; print(r'   omega    = %2.3e  rad/s ' % s)
s = 2*pi/w; print(r'   period T  = %2.3e  s ' % s)
s = L1*1e9; print(r'   deBroglie lambda1 = %2.3f nm ' % s)
s = real(L2*1e9); print(r'   deBroglie lambda2 = %2.3f nm '  % s)
print(r'   Jinc    = %2.2e  1/s ' % Jinc)
s = real(Jref);print(r'   Jref    = %2.2e 1/s ' % s)
s = real(Jinc+Jref); print(r'   Jnet    = %2.2e  1/s ' % s)
s = real(Jtrans); print(r'   Jtrans  = %2.2e  1/s ' % s)
if E0>U0:
   print(r'   Reflection coeff.   R = %2.2f  ' % Rc)
   print(r'   Transmission coeff. T = %2.2f  ' % Tc)


#%% GRAPHICS
plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (6,4)

fig1, axes = plt.subplots(nrows=2, ncols=1)

R = 0
axes[R].set_ylabel('$\psi$',color = 'black',fontsize = 14)

axes[R].xaxis.grid()
axes[R].yaxis.grid()
axes[R].set_xlim(xMin*1e9,xMax*1e9)
axes[R].set_xticks(np.arange(-0.6,0.7,0.2))
xP = x1*1e9; yP = psi1R
axes[R].plot(xP, yP,'b', lw = 2)
xP = x2*1e9; yP = psi2R
axes[R].plot(xP, yP,'b', lw = 2)
xP = x1*1e9; yP = psi1I
axes[R].plot(xP, yP,'r', lw = 1)
xP = x2*1e9; yP = psi2I
axes[R].plot(xP, yP,'r', lw = 1)

R = 1
axes[R].set_xlabel('x  [ nm ] ',color = 'black',fontsize = 10)
axes[R].set_ylabel('|$\psi|^2$',color = 'black',fontsize = 14)
axes[R].xaxis.grid()
axes[R].yaxis.grid()
axes[R].set_xlim(xMin*1e9,xMax*1e9)
xP = x1*1e9; yP = probD1
axes[R].plot(xP, yP,'b', lw = 2)
xP = x2*1e9; yP = probD2
axes[R].plot(xP, yP,'b', lw = 2)
axes[R].fill_between(x1*1e9, probD1,color = [1,0,1],alpha=0.2)
axes[R].fill_between(x2*1e9, probD2,color = [1,1,0],alpha=0.2)

fig1.tight_layout()

fig1.savefig('a1.png')


#%%
# Reflection and transmission coefficients
U = 50
N = 299
E = linspace(U,250,N)
k1 = sqrt(2*me*E)/hbar           # propagation constant x < 0    [1/m]
k2 = sqrt(2*me*(E-U))/hbar       # propagation constant x > 0    [1/m]
Tc = 4*k1*k2/(k1+k2)**2

Rc = 1 - Tc



plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (4,3)
fig, ax = plt.subplots(1)
#fig, plt.rcParams['font.size'] = 12
# fig, plt.rcParams["figure.figsize"] = (3,3)
ax.xaxis.grid()
ax.yaxis.grid()
ax.set_ylabel('R, T',color= 'black')
ax.set_xlabel('E/U',color = 'black')
ax.set_xlim([0, 5])
ax.set_ylim([0, 1.1])
#ax.set_xticks(np.arange(0,101,20))
#ax.set_yticks(np.arange(-20,81,20))

xP = E/U; yP = Tc
ax.plot(xP,yP,'b',lw = 2, label = 'T')
yP = Rc
ax.plot(xP,yP,'r',lw = 2, label = 'R')

xP = [0,1]; yP = [0,0]
ax.plot(xP,yP,'b',lw = 2)

xP = [0,1]; yP = [1,1]
ax.plot(xP,yP,'r',lw = 2)


ax.legend()
fig.tight_layout()

fig.savefig('a2.png')

