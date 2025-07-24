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
def lorenz(t, state): 
    x, y, z = state
    # I = 0
    # if t >= T1: I = I0
    # if t > T2: I = 0
    dx = y + 3*x**2 - x**3 - z + I
    dy = 1 - 5*x**2 - y
    dz = r*(s*(x-xc) - z)
    return [dx, dy, dz]  


#%% INPUTS >>>
# Initial conditions
x1,y1,z1 = 0.1, 1.0, 0.2

# Adpation current
I,s = 3, 4  
# Time span for solving ODE
t1 = 0; t2 = 1400; nT = 9999     
t = linspace(t1,t2,nT)         
# External current stimulus
nR = 1200
R = linspace(0.005,0.05,nR)    
# Initial conditions
xc = -0.5*(1+sqrt(5)) 
xc = -1.6
u0 = [xc,y1,z1]

#%% Graphics setup
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,4)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r',fontsize = 12)
ax.set_ylabel('ISI',fontsize = 12)
#ax.set_xlim([1,4]); ax.set_ylim([0,200])
ax.grid()


#%% Solve ODEs 
for c in range(nR):
    r = R[c]
    sol = odeint(lorenz, u0, t, tfirst=True)
    xS = sol[:,0]     

    q = find_peaks(xS, height = 1)
    q = q[0]
    tB = t[q]
    L1 = int(len(tB)/1.5); L2 = round(len(tB))
    tBF = tB[L1:L2]
    L = len(tBF)
    dt = zeros(L) 
    for q in range(L-1):
        dt[q] = tBF[q+1]-tBF[q]

    dt[-1] = dt[-2]

    Rp = np.ones(len(dt))*R[c]
    ax.plot(Rp,dt,'bo',ms = 0.3)
    
    fig1.tight_layout()



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
#