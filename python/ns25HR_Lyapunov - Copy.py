# -*- coding: utf-8 -*-
'''
AA\ns\pyCode\ns25_HR.py      July 2025
# COMPLEX SYSTEMS: NEUROSCIENCE
# HINDMARSH-ROSE MODEL BURSTING NEURONS

#Ian Cooper
# https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation
#    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/ns25_HR.pdf

https://pmc.ncbi.nlm.nih.gov/articles/PMC9871131/

https://scicomp.stackexchange.com/questions/36013/numerical-computation-of-lyapunov-exponent


'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, real, imag, sqrt 
from numpy.linalg import eig
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import random
from scipy.optimize import fsolve
from scipy.signal import find_peaks
from sympy import symbols, Eq, solve

tStart = time.time()

plt.close('all')


#%% FUNCTIONS  Solve ODE for x,y    
# def lorenz(t, state): 
#     x, y, z = state
#     dx = y + 3*x**2 - x**3 - z + I
#     dy = 1 - 5*x**2 - y
#     dz = r*(s*(x-xc) - z)
#     return [dx, dy, dz]  


# #%% INPUTS >>>
# # Initial conditions
# x1,y1,z1 = 0.1, 1.0, 0.2
# xc = -0.5*(1+sqrt(5))
# # Adpation current
# r,s = 0.005, 4  
# # Time span for solving ODE
# t1 = 0; t2 = 2000; nT = 9999     
# t = linspace(t1,t2,nT)         
# # External current stimulus
# nI = 599
# Iext = linspace(1,4,nI)    

# #%% Graphics setup
# plt.rcParams['font.size'] = 12
# plt.rcParams["figure.figsize"] = (6,4)
# fig1, ax = plt.subplots(nrows=1, ncols=1)
# ax.set_xlabel('LE',fontsize = 12)
# ax.set_ylabel('ISI',fontsize = 12)
# ax.set_xlim([1,4]); #ax.set_ylim([0,200])
# ax.grid()



# #%% Solve ODEs 
# L = zeros(nI)
# for c in range(nI):
#     I = Iext[c]
#     x1,y1,z1 = 0.1, 1.0, 0.2
#     u0 = [x1,y1,z1]
#     sol = odeint(lorenz, u0, t, tfirst=True)
#     xS1 = sol[:,0] 
#     x2 = x1 + 0.0001    
#     u0 = [x2,y1,z1]
#     sol = odeint(lorenz, u0, t, tfirst=True)
#     xS2 = sol[:,0]
#     S = sqrt((xS1-xS2)**2)
#     S1 = S[6001:9990]
#     S2 = S[6000:9989]
#     S3 = np.log(S1/S2)
#     S4 = sum(S3)
#     L[c] = S4
#     ax.plot(I,S4,'bo',ms = 0.5)


#%%
plt.close('all')
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,4)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('real(x$_C$)',fontsize = 12)
ax.set_ylabel('imag(x_C$)',fontsize = 12)
ax.set_xlim([-2,2]); ax.set_ylim([-2,2])
ax.grid()
fig1.tight_layout()

S = 0.6
r,s = 0.002, 4
xR = xc = -0.5*(1+sqrt(5)) 
Imax = 20
I0 = linspace(0,Imax,299)

for c in range(299):
   c1 = 1; c2 = 2; c3 = r*s*xR; c4 = + r*s*xR - 1 - I0[c]
   coeff = [c1,c2,c3,c4]
   xC = np.roots(coeff)
   yC = 1 - 5*xC**2  
   zC = r*s*xC - r*s*xR 
   ax.plot(real(xC[0]),imag(xC[0]),'bo',ms = S)
   ax.plot(real(xC[1]),imag(xC[1]),'ro',ms = S)
   ax.plot(real(xC[2]),imag(xC[2]),'mo',ms = S)
   if c == 0:
      ax.plot(real(xC[0]),imag(xC[0]),'bo',ms = 8)
      ax.plot(real(xC[1]),imag(xC[1]),'ro',ms = 8)
      ax.plot(real(xC[2]),imag(xC[2]),'mo',ms = 8)   
#print(xC,yC,zC)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,4)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('I',fontsize = 12)
ax.set_ylabel('x$_C$',fontsize = 12)
ax.set_xlim([0,Imax]); ax.set_ylim([-2,2])
ax.grid()
fig2.tight_layout()

S = 2
for c in range(299):
   c1 = 1; c2 = 2; c3 = r*s*xR; c4 = + r*s*xR - 1 - I0[c]
   coeff = [c1,c2,c3,c4]
   xC = np.roots(coeff)
   yC = 1 - 5*xC**2  
   zC = r*s*xC - r*s*xR 
   ax.plot(I0[c], real(xC[0]),'bo',ms = S)
   ax.plot(I0[c], imag(xC[0]),'k+',ms = S)
   #ax.plot(I0[c], real(xC[1]),'ro',ms = S)
   #ax.plot(I0[c], imag(xC[1]),'m+',ms = S)
   #ax.plot(I0[c], real(xC[2]),'bo',ms = S)
   #ax.plot(I0[c], imag(xC[2]),'r+',ms = S)
   
   #ax.plot(real(xC[1]),imag(xC[1]),'ro',ms = S)
   #ax.plot(real(xC[2]),imag(xC[2]),'mo',ms = S)



#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)

#%%
'''
fig1.savefig('a1.png')
fig2.savefig('a2.png')
fig3.savefig('a3.png')
fig4.savefig('a4.png')
fig5.savefig('a5.png')



'''
