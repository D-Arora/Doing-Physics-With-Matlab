# cs_003.py
# Ian Cooper
# Jan 2024
# COMPLEX SYSTEMS
#   DYNAMICAL SYSTEMS WITH TWO DEGREES OF FREEDOM:
#    PREDATOR-PREY SYSTEMS (Lotka-Volterra Equations)

# # Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_003.pdf

# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.signal import chirp, find_peaks, peak_widths
import time
#%%

# SETUP
N = 9999
R = np.zeros(N)       # rabbits
W = np.zeros(N)       # wolfs
R[0] = 1200; W[0] = 50
a = np.zeros(4)       # constants
a[0] = 8e-2; a[1] = 1e-3; a[2] = 2e-2; a[3] = 2e-5

# Time span   [months]
tMin = 0; tMax = 500
t = np.linspace(tMin,tMax,N)
h = t[2] - t[1]

# Solve predator-prey equations
for c in range(N-1):
    R[c+1] = R[c] + h*( a[0]*R[c] - a[1]*R[c]*W[c] )
    W[c+1] = W[c] + h*( -a[2]*W[c] + a[3]*R[c+1]*W[c] )

# Equilibrium: steat-sate values
Rss = a[2]/a[3]
Wss = a[0]/a[1]

# Peaks
peaks, _ = find_peaks(R)
TR = t[peaks]
pR = 0.5*(TR[2]-TR[0])
peaks, _ = find_peaks(W)
TW = t[peaks]
pW = 0.5*(TW[2]-TW[0])
RW = (TW[0] - TR[0])/pR
peakRW = (TW[0] - TR[0])/pR
# Console output
print(' ')
print('rabbits  Rss = %2.0f' % Rss  )
print('wolfs    Wss = %2.0f' % Wss  )
print('rabbits  Rmax = %2.0f' % max(R), '  Rmin = %2.0f' % min(R)  )
print('wolfs    Wmax = %2.0f' % max(W), '  Wmin = %2.0f' % min(W) )
print(' ')
print('rabbits   period =  %2.2f  months' % pR  )
print('wolfs     period = %2.2f   months'  % pW  )
print('peaks: phase R --> W = %2.2f  cycles'  % peakRW  )

#%%
#GRAPHICS  ===============================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = ['Garamond']
plt.rcParams['font.size'] = 14
font1 = {'color':'black','size':14} 


# Figure 1     R vs t
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.2,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  
col = [0,0,1]

xP = t
yP = R
plt.plot(xP,yP,linewidth=2,color = col)

yR = np.arange(0,3001,500) 
plt.yticks(yR)
xR = np.arange(0,501,100) 
plt.xticks(xR)

plt.grid('visible')
plt.xlabel("t  [ months ]", fontdict = font1)
plt.ylabel('rabbits  R ', fontdict = font1)

plt.savefig('a_cs_001.png')    

# Figure 2      W vs t
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.2,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  
col = [1,0,0]

xP = t
yP = W
plt.plot(xP,yP,linewidth=2,color = col)
xR = np.arange(0,501,100)
plt.xticks(xR) 
yR = np.arange(0,151,25) 
plt.yticks(yR)

plt.grid('visible')
plt.xlabel("t  [ months ]", fontdict = font1)
plt.ylabel('wolfs   W ', fontdict = font1)
plt.savefig('a_cs_002.png')

# Figure 3     phase space diagram

fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.2,\
                    right = 0.92, hspace = 0.60,wspace=0.5)  
col = [0,0,1]

xP = R
yP = W
plt.plot(xP,yP,linewidth=2,color = col)

plt.plot(Rss,Wss,'or',ms = 8)
plt.plot(R[0],W[0],'og',ms = 8)

yR = np.arange(0,141,20) 
plt.yticks(yR)
xR = np.arange(0,3001,500) 
plt.xticks(xR)

plt.grid('visible')
plt.ylabel("wolfs  W ", fontdict = font1)
plt.xlabel('rabbits  R ', fontdict = font1)

plt.savefig('a_cs_003.png')  

# figure 4   Quiver plot 
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.2,\
                    right = 0.92, hspace = 0.60,wspace=0.5)  
col = [0,0,1]

X = np.linspace(100, 2900, 9)
Y = np.arange(10,140,9)
 
XX, YY = np.meshgrid(X, Y)
u = a[0]*XX - a[1]*XX*YY
v = -a[2]*YY + a[3]*XX*YY

us = u/(u**2 + v**2)**0.5
vs = v/(u**2 + v**2)**0.5
plt.quiver(XX, YY, us, vs,color = [col])
plt.plot(Rss,Wss,'or',ms = 8)

xR = np.arange(0,3001,500)
plt.xticks(xR) 
yR = np.arange(0,141,20) 
plt.yticks(yR)

plt.ylabel("wolfs  W ", fontdict = font1)
plt.xlabel('rabbits  R ', fontdict = font1)

plt.savefig('a_cs_004.png')

# figure 5   Stream plot 
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.2,\
                    right = 0.92, hspace = 0.60,wspace=0.5)  
col = [0,0,1]

X = np.linspace(100, 2900, 9)
Y = np.arange(10,140,9)
 
XX, YY = np.meshgrid(X, Y)
u = a[0]*XX - a[1]*XX*YY
v = -a[2]*YY + a[3]*XX*YY

us = u/(u**2 + v**2)**0.5
vs = v/(u**2 + v**2)**0.5
plt.streamplot(XX, YY, us, vs)
plt.plot(Rss,Wss,'or',ms = 8)

xR = np.arange(0,3001,500)
plt.xticks(xR) 
yR = np.arange(0,141,20) 
plt.yticks(yR)

plt.ylabel("wolfs  W ", fontdict = font1)
plt.xlabel('rabbits  R ', fontdict = font1)

plt.savefig('a_cs_005.png')


#%%  figure 6    R and W vs t

fig, ax = plt.subplots(figsize = (4, 3))

ax2 = ax.twinx()
ax.plot(t, R, color = 'b')
ax2.plot(t, W, color = 'r')

xR = np.arange(0,501,100) 
ax.set_xticks(xR)
plt.grid('visible')
yR = np.arange(0,3001,500)
ax.set_yticks(yR)
yR = np.arange(0,151,25)
ax2.set_yticks(yR)

ax.grid(axis = "x")
ax.tick_params(axis='y', colors='blue') 
ax2.tick_params(axis='y', colors='red') 

ax.set_xlabel('t  [ months ]') 
ax.set_ylabel('rabbits  R', color = 'b')
ax2.set_ylabel('wolfs  W', color = 'r')
  
plt.tight_layout()      # defining display layout 
plt.show()

plt.savefig('a_cs_006.png')


#%%
# Figure 7     quiver plot and trajectories

# figure 4   Quiver plot 
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.2,\
                    right = 0.92, hspace = 0.60,wspace=0.5)  


col = [1,0,1]
xP = R
yP = W
plt.plot(xP,yP,linewidth=2,color = col)
plt.plot(Rss,Wss,'or',ms = 6)
plt.plot(R[0],W[0],'og',ms = 6)    
    
col = [0,0,1]
X = np.linspace(100, 2900, 9)
Y = np.arange(10,140,9)
 
XX, YY = np.meshgrid(X, Y)
u = a[0]*XX - a[1]*XX*YY
v = -a[2]*YY + a[3]*XX*YY

us = u/(u**2 + v**2)**0.5
vs = v/(u**2 + v**2)**0.5
plt.quiver(XX, YY, us, vs,color = [col])
plt.plot(Rss,Wss,'or',ms = 8)

xR = np.arange(0,3001,500)
plt.xticks(xR) 
yR = np.arange(0,141,20) 
plt.yticks(yR)

plt.ylabel("wolfs  W ", fontdict = font1)
plt.xlabel('rabbits  R ', fontdict = font1)

plt.savefig('a_cs_007.png')

