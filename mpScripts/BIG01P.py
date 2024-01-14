# BIG01P.py
# Ian Cooper
# Jan 2024
# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/BIG-A.pdf
# Beta-cell mass B, insulin I, glucose G dynamics
# Fixed points: Bf If Gf
# B-cell mass dynamics: G(B) I(B) B plots and fixed points
# (dB/dt)/B vs G plot

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
#from pylab import *
from math import *
from numpy import linspace,matrix,array
from numpy import *
from scipy.linalg import *
from scipy.optimize import fsolve


#%%
# CONSTANTS  ==========================================================
k = np.zeros(13)
k[1]  = 864         # 864
k[2]  = 1.44        # 1.44
k[3]  = 1        # 1
k[4]  = 0.01        # 0.01
k[5]  = 40         # 40
k[6]  = 20000       # 20000
k[7]  = 432         # 432
k[8]  = 0.06        # 0.06
k[9]  = 10e-4       # 10e-4
k[10] = 2.4e-6      # 2.4e-6
k[11] = 3.49        # 3.49
k[12] = np.pi/3050   # pi/3050



#%%
# FIXED POINTS: Steady values Gss, Iss, Bss as a function of k1
# Must select starting values carefully to obtain correct solution

def myFunction(u):
   x = u[0]
   y = u[1]
   z = u[2]
   F = np.empty((3))
   F[0] = p - k[2]*x - k[3]*y*x/(k[4]*x+1)
   F[1] = z*k[5]*x**2 / (k[6]+x**2)-k[7]*y
   F[2] = (-k[8] + k[9]*x - k[10]*x**2)*z - k[11]*np.tan(k[12]*z)
   return F


# >>>  Starting values Gss, Iss, Bss
zGuess = np.array([94,38,1500]) 

nP = 999
Gss = np.zeros(nP); Iss = np.zeros(nP); Bss = np.zeros(nP)
P = linspace(200,5200,nP)     # k1

for c in range(nP):
    p = P[c]
    z = fsolve(myFunction,zGuess)
    Gss[c] = z[0]
    Iss[c] = z[1]
    Bss[c] = z[2]
    ZGuess = np.array([Gss[c],Iss[c],Bss[c]])
    

#%%  GRAPHICS  

def plotFn(xP,yP,limX,tX,limY,tY,txtx,txtY):
    font1 = {'family':'Cambria','color':'black','size':12}
    plt.rcParams['font.family'] = ['Tahoma']
    plt.rcParams['font.size'] = 12

    fig = plt.figure(figsize=(5.5,3))
    fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.15,\
                    right = 0.94, hspace = 0.45,wspace=0.5)
     
    plt.plot(xP,yP,linewidth=2,color='b')
    plt.xlim(limX)
    plt.xticks(tX)
    plt.ylim(limY)
    plt.yticks(tY)
    plt.grid('visible')  
    plt.xlabel(txtX,fontdict = font1)
    plt.ylabel(txtY, fontdict = font1)
    return
    
#%%  Fig 1:  Gss  vs k1
y1 = 0; y2 = 200; dy = 25
x1 = 000; x2 = 5150; dx = 500
limX = [x1, x2]
tX = range(x1,x2,dx)
limY= [y1, y2]
tY = range(y1,y2,dy)
txtX = '$k_{1}$    $[ mg.dL^{-1}.d^{-1} ]$ '
txtY = '$G_{ss}$    $[ mg.dL^{-1} ]$ '
xP = P
yP = Gss
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)

#%%   Fig 2: Iss  vs k1
y1 = 0; y2 = 85; dy = 20
x1 = 000; x2 = 5150; dx = 500
limX = [x1, x2]
tX = range(x1,x2,dx)
limY= [y1, y2]
tY = range(y1,y2,dy)
txtX = '$k_{1}$    $[ mg.dL^{-1}.d^{-1} ]$ '
txtY = '$I_{ss}$    $[ mU.L^{-1} ]$ '
xP = P
yP = Iss
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)

#%%  Fig 3: Bss vs k1s
y1 = 0; y2 = 1550; dy = 250
x1 = 000; x2 = 5150; dx = 500
limX = [x1, x2]
tX = range(x1,x2,dx)
limY= [y1, y2]
tY = range(y1,y2,dy)
txtX = '$k_{1}$    $[ mg.dL^{-1}.d^{-1} ]$ '
txtY = '$B_{ss}$    $[ mg ]$ '
xP = P
yP = Bss
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)

#%% Fig 4:  Michaelis function   GH vs G
GM = linspace(0,600,299)
M = GM/(k[4]*GM+1)
y1 = 0; y2 = 105; dy = 20
x1 = 0; x2 = 605; dx = 100
limX = [x1, x2]
tX = np.arange(x1,x2,dx)
limY= [y1, y2]
tY = np.arange(y1,y2,dy)
txtX = '$G$    $[ mg.dL^{-1} ]$ '
txtY = '$G_M}$ $[mg.dL^{-1}]$ '
xP = GM
yP = M
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)

#%% Fig 5:  Hill function
GH = linspace(0,600,299)
k6 = 150
H = GH**2/(k6**2 + GH**2)
y1 = 0; y2 = 1.05; dy = 0.2
x1 = 0; x2 = 605; dx = 100
limX = [x1, x2]
tX = np.arange(x1,x2,dx)
limY= [y1, y2]
tY = np.arange(y1,y2,dy)
txtX = '$G$    $[ mg.dL^{-1} ]$ '
txtY = '$G_H}$ '
xP = GH
yP = H
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)
