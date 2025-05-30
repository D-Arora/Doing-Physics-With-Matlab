# BIG013.py
# Ian Cooper
# Jan 2024
# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/BIG-A.pdf
# Beta-cell mass B, insulin I, glucose G dynamics
# Gss Iss Bss dG/dt  vs k5    insulin secretion rate

# LIBRARIES  ==============================================================
import numpy
import numpy as np
import scipy.special 
import scipy.special as sps
from scipy import *
import scipy.linalg as la
import math
import os
import matplotlib.pyplot as plt
from pylab import *
from math import *
from numpy import linspace,matrix,array
from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve
from scipy.integrate import odeint, solve_ivp
#o s.system('cls')

# FUNCTIONS  ==============================================================
#beta-cell mass (z) insulin (y)   glucose (x)  
def lorenz(t, state, k):    
    x, y, z = state
    
    dx = k[1] - k[2]*x - k[3]*x*y/(k[4]*x+1) 
    dy = k[5]*z*x**2/(k[6]+x**2) - k[7]*y
    dz = (-k[8] + k[9]*x - k[10]*x**2 )*z - k[11]*np.tan(k[12]*z)
        
    return [dx, dy, dz]



 #%%   
# CONSTANTS  ===============================================================
k = np.zeros(13)
k[1]  = 864    # 864
k[2]  = 1.44     # 1.44
k[3]  = 1  # 1
k[4]  = 0.01     # 0.01
k[5]  = 40 # 40
k[6]  = 20000    # 20000
k[7]  = 432      # 432
k[8]  = 0.06    # 0.06
k[9]  = 10e-4  # 10e-4
k[10] = 2.4e-6   # 2.4e-6
k[11] = 3.49      # 3.49
k[12] = np.pi/3050    # pi/3050


#%%
# SOLVE ODEs  ============================================================ 
y0 = [80, 17, 755]  # Inital Conditions: G0  I0   B0
N = 999            # Time interval for t
t1 = 0
t2 = 500
tSpan = np.linspace(t1,t2,N)
t = tSpan

wN = 599
Gss = np.zeros(wN)
Iss = np.zeros(wN)
Bss = np.zeros(wN)
Gs = np.zeros(wN)
Gp = np.zeros(wN)
w = linspace(3.7,50,wN)

for n in range(wN):
    k[5] = w[n]
    sol = odeint(lorenz, y0, tSpan, args = (k,), tfirst=True)
    G = sol[:,0]
    I = sol[:,1]
    B = sol[:,2]
    Gss[n] = G[N-1]
    Iss[n] = I[N-1]
    Bss[n] = B[N-1]        

    B5 = Bss[n]
    a = k[10]*B5
    b = -k[9]*B5
    c = k[8]*B5 + k[11]*np.tan([k[12]*B5])

    Gs[n] = (-b + (b**2 - 4*a*c)**0.5)/(2*a)
    Gp[n] = (-b - (b**2 - 4*a*c)**0.5)/(2*a)
    
    gradGp = abs(np.gradient(Gp))
    
    
# GRAPHICS  ================================================================

fig1 = plt.figure(figsize=(4,3))
fig1.subplots_adjust(top=0.92, bottom = 0.20, left = 0.2,\
                    right = 0.96, hspace = 0.45,wspace=0.5)
xP = w
yP = Gs
plt.plot(xP,yP,'r',linewidth = 2)
yP = Gp
plt.plot(xP, yP,'b',linewidth = 2)
xP = [0, 50]
yP = [100, 100]
plt.plot(xP,yP,'m')
yP = [125,125]
plt.plot(xP,yP,'m')   
plt.grid('visible')
plt.ylabel('G$_{ss}$ ',fontsize = 14)
plt.xlabel('k$_{5}$',fontsize = 14)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.text(2,70,'NGT')
plt.text(42,105.5,'IGT')
plt.text(42,130,'T2D')
plt.show()

fig2 = plt.figure(figsize=(4,3))
fig2.subplots_adjust(top=0.92, bottom = 0.20, left = 0.2,\
                    right = 0.96, hspace = 0.45,wspace=0.5)
xP = w
yP = gradGp
plt.plot(xP,yP,'b',linewidth = 2)
plt.grid('visible')
plt.ylabel('dG$_{ss}$ / dk$_{5}$ ',fontsize = 14)
plt.xlabel('k$_{5}$',fontsize = 14)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.show()

fig3 = plt.figure(figsize=(4,3))
fig3.subplots_adjust(top=0.92, bottom = 0.20, left = 0.2,\
                    right = 0.96, hspace = 0.45,wspace=0.5)
xP = w
yP = Iss
plt.plot(xP,yP,'b',linewidth = 2)
plt.grid('visible')
plt.ylabel('I$_{ss}$ ',fontsize = 14)
plt.xlabel('k$_{5}$',fontsize = 14)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.show()


fig4 = plt.figure(figsize=(4,3))
fig4.subplots_adjust(top=0.92, bottom = 0.20, left = 0.2,\
                    right = 0.96, hspace = 0.45,wspace=0.5)

xP = w
yP = Bss
plt.plot(xP,yP,'b',linewidth = 2)
plt.grid('visible')
plt.ylabel('B$_{ss}$ ',fontsize = 14)
plt.xlabel('k$_{5}$',fontsize = 14)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.show()