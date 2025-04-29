# -*- coding: utf-8 -*-
"""

April 2025

emRSBessel02.py               April 2025
 

COMPUTATIONAL OPTICS
   This is an example Code for computing the Hermite-Gaussian function
   w0 = 1    2x the standard deviation of the gaussian

Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/
 Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/emRS02.pdf

https://dcc.ligo.org/public/0091/G1200548/001/G1200548-v1.pdf

"""

import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, sqrt 
import matplotlib.pyplot as plt
import time
from scipy import special

tStart = time.time()

#%%
n = 3
H = special.hermite(n, monic=False)
print('Hermite polynomial,coefficients, H(1)')
print(H)
print(H(1))
x = np.linspace(-3, 3, 400)
y = H(x)

plt.rcParams['font.size'] = 10
plt.rcParams["figure.figsize"] = (5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(x, y,'b',lw = 2)
ax.set_title("Hermite polynomial of degree n = %0.0f" %n,fontsize = 10)
ax.set_xlabel("x")
ax.set_ylabel("H$_n$(x)")
ax.grid()
fig1.tight_layout()


#%%  HERMITRE-GAUSSS [2D] FUNCTION
m = 3           # x and y orders
n = 3    
w0 = 1          # standard eviation of Gausssian function
q = sqrt(2)/w0  # constant
# XY Grid 
x = np.linspace(-3, 3, 199)
y = x
X,Y = np.meshgrid(x,y)    
# Gaussian function
E0 = exp(-(X**2 + Y**2)/w0**2)
# XY Hermite functions
xH = special.hermite(m, monic=False)
yH = special.hermite(n, monic=False)
Hx  = xH(q*X); Hy = yH(q*Y)
# XY Hermite-Gaussian function
HH = Hx*Hy*exp(-X**2/w0)*exp(-Y**2/w0)
HH = HH/amax(HH)

#%% GRAPHICS
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (4.2,4.2)

fig2, ax = plt.subplots(nrows=1, ncols=1)
# fig2.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.0,\
#                     right = 1, hspace = 0.10,wspace=0.1)
    
cf = ax.pcolor(X,Y,HH, cmap='jet')      # Greys_r

ax.set_aspect('equal', adjustable='box')
ax.set_title('m = %0.0f   '%m + 'n = %0.0f' %n,fontsize = 8)
ax.set_xticks([]); ax.set_yticks([])
fig2.colorbar(cf, ax=ax)
fig2.tight_layout()

plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (4.2,4.2)

fig3, ax = plt.subplots(nrows=1, ncols=1)
# fig2.subplots_adjust(top = 0.92, bottom = 0.23, left = 0.0,\
#                     right = 1, hspace = 0.10,wspace=0.1)
    
cf = ax.pcolor(X,Y,HH**2, cmap='jet')      # Greys_r

ax.set_aspect('equal', adjustable='box')
ax.set_title('m = %0.0f   '%m + 'n = %0.0f' %n,fontsize = 8)
ax.set_xticks([]); ax.set_yticks([])
fig3.colorbar(cf, ax=ax)
fig3.tight_layout()


#%%
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')

