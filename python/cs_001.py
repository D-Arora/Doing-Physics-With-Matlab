# cs_001.py
# Ian Cooper
# Jan 2024
# COMPLEX SYSTEMS
#   Simulating discrete-time models with one variable
# # Website: https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_001.pdf


# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

plt.close('all')

#%% CELL 1   Exponetial growth
# SETUP ====================================================================
a = 1.2            # <<<<<  constant
x0 = 1            # <<<<<  intial condition
nT = 21          # <<<<<  number of time steps

x = np.zeros(nT)     # x variable
x[0] = x0
R = np.arange(nT)

# Iteriation
for c in range(nT-1):
    x[c+1] = a*x[c]
    
b = np.log(x[nT-1]/x0)/(nT-1)    
y = x0*np.exp(b*R)  

    
# GRAPHICS  ===============================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = ['Garamond']
plt.rcParams['font.size'] = 14
font1 = {'color':'black','size':14}

fig = plt.figure(figsize = (5, 3))
fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.15,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  

# Figure 1
xP = R
yP = x
plt.plot(xP,yP,linewidth=2,color = 'b')
yP = y
plt.plot(xP,yP,'o',color = 'r',ms = 4)
yR = np.arange(0,41,5) 
plt.yticks(yR)
xR = np.arange(0,21,5) 
plt.xticks(xR)
plt.grid('visible')
plt.ylabel("x(c+1)", fontdict = font1)
plt.xlabel('x(c)', fontdict = font1)
plt.savefig('a_cs_001.png')    

#%% CELL 2    Fish population with harvesting
# SETUP ====================================================================
a = 1.1         # <<<<<  constant          a > 0
b = 9       # <<<<< harvesting constant   b > 0
x0 = 100           # <<<<<  intial condition
nT = 250        # <<<<<  number of time steps

col = [0,0,0]    # <<<<< Color line in plot

x = np.zeros(nT)     # x variable
x[0] = x0
R = np.arange(nT)

# Iteriation
for c in range(nT-1):
    x[c+1] = a*x[c] - b
    
x[x<0] = 0


#%%    
# GRAPHICS  ===============================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = ['Garamond']
plt.rcParams['font.size'] = 14
font1 = {'color':'black','size':14}

fig = plt.figure(figsize = (5, 3))
fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.15,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  

# Figure 2
xP = R
yP = x
plt.plot(xP,yP,linewidth=2,color = col)
yR = np.arange(0,501,100) 
#plt.yticks(yR)
xR = np.arange(0,21,5) 
#plt.xticks(xR)
plt.grid('visible')
plt.ylabel('fish population', fontdict = font1)
plt.xlabel('t', fontdict = font1)
plt.savefig('a_cs_002.png')    
