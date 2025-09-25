# cs_002.py
# Ian Cooper
# Jan 2024
# COMPLEX SYSTEMS
#   Simulating discrete-time models with one variable
# # Website: https://d-arora.github.io/Doing-Physics-With-Matlab/
# Documentation: https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cs_001.pdf
# Newton's law of cooling
#  CSI:   MURDER - TIME OF DEATH ?

# LIBRARIES  ================================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#%%

def temp(T0,Tenv,R,N):
    for c in range(N-1):
        T[c+1] = T[c] - h*R*( T[c]-Tenv ) 
    return T

# >>> INPUTS
# Line color for plot
col = [0,0,1]
# Environmental temperature
Tenv = 21
# R found by a trail-and-error process  Tend = 28.5 in 60 minutes
R = 6.38e-3


# SETUP
# 
T0 = 37     # Initial body temperature
T1 = 32     # Body temperature when dead body found at 5:00 pm
T2 = 28.5   # Body temperature 1 hour after found  6:00 pm

# SETUP
# Time span
tMin = 0; tMax = 140; N = 999
t = np.linspace(tMin,tMax,N)
h = t[2] - t[1]

# Temperature
T = np.zeros(N)
T[0] = T0
T = temp(T1,Tenv,R,N) 

# Target time interval:  eant dt = 60
t1 =  t[T<T1][0]
t2 = t[T<T2][0]
dt =  t2 - t1  


#%%
# GRAPHICS  ===============================================================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = ['Garamond']
plt.rcParams['font.size'] = 14
font1 = {'color':'black','size':14}

fig = plt.figure(figsize = (5, 4))
fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.15,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  

# Figure 1
xP = t
yP = T
plt.plot(xP,yP,linewidth=2,color = col)

xP = [0, tMax]; yP = [T1, T1]
plt.plot(xP,yP,linewidth=1,color = [0,0,0])
yP = [T2, T2]
plt.plot(xP,yP,linewidth=1,color = [0,0,0])

xP = [t1, t1]; yP = [25, 40]
plt.plot(xP,yP,linewidth=1,color = [0,0,0])
xP = [t2, t2]
plt.plot(xP,yP,linewidth=1,color = [0,0,0])

plt.ylim([26,40])
yR = np.arange(26,41,2) 
plt.yticks(yR)
plt.xlim([0,tMax])
xR = np.arange(0,tMax+12,20) 
plt.xticks(xR)

plt.grid('visible')
plt.xlabel("t  [ min ]", fontdict = font1)
plt.ylabel('T  [ deg C ]', fontdict = font1)


plt.text(87,27,'Tenv = %2.0f'  % Tenv)
plt.text(29.2,39,'t1 = %2.1f'  % t1)
plt.text(89.2,39,'t2 = %2.1f'  % t2)
plt.text(81,33,'dt = %2.2f'    % dt)
#plt.title("Target Electron Temperature =%1.0f" %Te_t + "[ev] \nTarget Density=%1.1f"%n_t + "[m^-3]")
plt.title("R = %2.3e" % R, fontsize = 12)
#plt.title("dt1 = {dt1} ".format(dt1 = dt1))

#lt.title("Target Electron Temperature={Te_t}[ev] \nTarget Density={n_t},[m^-3]".format(Te_t=Te_t, n_t=n_t))

plt.savefig('a_cs_001.png')    


#%%
# Figure 2
fig = plt.figure(figsize = (5, 4))
fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.15,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  

xP = t
yP = T
plt.plot(xP,yP,linewidth=2,color = col)

plt.ylim([26,40])
yR = np.arange(26,41,2) 
plt.yticks(yR)
plt.xlim([0,tMax])
xR = np.arange(0,tMax+12,20) 
plt.xticks(xR)

plt.grid('visible')
plt.xlabel("t  [ min ]", fontdict = font1)
plt.ylabel('T  [ deg C ]', fontdict = font1)


plt.savefig('a_cs_002.png')    



# Estimate# d R values from trial-and-error process

# Tenv = 21  R = 6.38e-3
# Tenv = 22  R = 7.17e-3
# Tenv = 23  R = 8.20e-3
