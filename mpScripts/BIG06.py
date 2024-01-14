# BIG06.py
# Ian Cooper
# Jan 2024
# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/BIG-A.pdf
# Beta-cell mass B, insulin I, glucose G dynamics
# Gss vs k10
# Population reduction term


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

# CONSTANTS  ===========================================================
p = np.zeros(12)
p[8] = 0.06               # k8
p[9] = 8.4e-4             # k9
p[10] = 2.4e-6            # k10

pMax = p[9]**2/(4*p[8])
p10 = arange(1e-6,pMax,0.01e-6)
L = len(p10)
Gfp = zeros(len(p10))
Gfs = zeros(len(p10))

# CALCULATIONS  ============================================================
# Steady-state glucse as a functin of k10
for c in range(L):
    Gfp[c] =  (p[9] - (p[9]**2-4*p10[c]*p[8])**0.5)/(2*p10[c])
    Gfs[c] =  (p[9] + (p[9]**2-4*p10[c]*p[8])**0.5)/(2*p10[c])


# GRAPHICS  ===============================================================
font1 = {'family':'Cambria','color':'black','size':14}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12
plt.rcParams['mathtext.default'] 

# Fig 1    Gss vs k10
fig1 = plt.figure(figsize=(5,3))
fig1.subplots_adjust(top=0.95, bottom = 0.29, left = 0.20,\
                    right = 0.94, hspace = 0.45,wspace=0.5)
xP = p10; yP = Gfp
plt.plot(xP,yP,linewidth=2,color='b',label = 'physiological')
xP = p10; yP = Gfs
plt.plot(xP,yP,linewidth=2,color='r',label = 'saddle')
#plt.xlim([0,40])
#plt.ylim([0,300])
#plt.xticks(np.arange(0,25,2))
plt.grid('visible')   
plt.xlabel('k$_{10}$   $[ mg^{2}.dL^{2}.d^{-1} ] $ ' ,fontdict = font1)
plt.ylabel('G$_{ss}$   [ mg.dL$^{-1}$ ]', fontdict = font1)
plt.title('  ')
plt.legend(fontsize = 10)
plt.show()



#%%   beta cell dynamics  reduction factor
b = np.linspace(0,1500,299)
k11 = 3.49
k12 = np.pi/3050

bR = k11*np.tan(k12*b)

fig2 = plt.figure(figsize=(4,3))
fig2.subplots_adjust(top=0.96, bottom = 0.20, left = 0.20,\
                    right = 0.94, hspace = 0.45,wspace=0.5)
font1 = {'family':'Cambria','color':'black','size':14}
    
xP = b; yP = bR
plt.plot(xP,yP,linewidth=2,color='b')

#plt.xlim([0,40])
#plt.ylim([0,300])
#plt.xticks(np.arange(0,25,2))
plt.grid('visible')   
plt.xlabel('$\\beta$' ,fontdict = font1)
plt.ylabel('$\\beta_{R}$ ', fontdict = font1)
plt.show()

 
