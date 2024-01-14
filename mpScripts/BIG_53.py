# -*- coding: utf-8 -*-

# BIG_51.py
# -*- coding: utf-8 -*-
# Ian Cooper
# 240107

# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Beta-cell mass B, insulin I, glucose G dynamics
# Time evolution of beta-cell mass, insulin, gluocose
# Solution of 3 coupled ODEs using the odeint function

# Response of glucose regulatory system after the consumption of food
# G-I index

# LIBRARIES  ==============================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
#from scipy.signal import find_peaks
from scipy.signal import chirp, find_peaks, peak_widths

#%%
# CONSTANTS  ==========================================================
k = np.zeros(13)
k[1]  = 864     # 864
k[2]  = 1.44       # 1.44
k[3]  = 1# 1
k[4]  = 0.01        # 0.01
k[5]  = 40          # 40
k[6]  = 150       # 20000
k[7]  = 432         # 432
k[8]  = 0.06        # 0.06
k[9]  = 10e-4       # 10e-4
k[10] = 2.4e-6      # 2.4e-6
k[11] = 3.49        # 3.49
k[12] = np.pi/3050   # pi/3050

col = [0,0,1]


# FIXED POINTS: Steady values Gss, Iss, Bss 
# Must select starting values carefully to obtain correct solution

def myFunction(u):
   x = u[0]
   y = u[1]
   z = u[2]
   F = np.empty((3))
   F[0] = k[1] - k[2]*x - k[3]*y*x/(k[4]*x+1)
   F[1] = z*k[5]*x**2 / (k[6]**2+x**2)-k[7]*y
   F[2] = (-k[8] + k[9]*x - k[10]*x**2)*z - k[11]*np.tan(k[12]*z)
   return F

zGuess = np.array([150,36,1400])     # Statring values Gss, Iss, Bss
z = fsolve(myFunction,zGuess)
Gss = z[0]
Iss = z[1]
Bss  = z[2]
Gph = k[1]/k[2]

print('Pathological fixed-point')
print('Bss = 0   Iss = 0    Gss = %5.2f' % Gph  )
print('  ')
print('Physiological fixed-point')
print('Bss = %5.2f ' % Bss,  'Iss = %5.2f' % Iss, 'Gss = %5.2f' % Gss  )
print('  ')


# SETUP  =================================================================
# Time span
N = 9999
t1 = 0
t2 = 1
t = np.linspace(t1,t2,N)
dt = t[2] - t[1]

# MEALS (STOMACH) ======================================================
# Inputs for 8 meals: Quantity G0 / Half life h > 0 / Meal time T
# >>> Meal quantity: glucose
G0 = np.array([100,200,300,00,00,00,0,0])  

# >>> Half-life [min] h > 0   
h  = np.array([20,20,20,20,20,20,1,1])  
tH = h/(60*24)
# >>> Time of meal  [hours]    
T  = np.array([2,2.25,3.0,4,8,16,0,0])   

M = 8    # max number of meals

kS = -np.log(1/2)/tH
T = T/24
Tind = T/dt
GSTO = np.zeros(N)
Gsto = np.zeros((N,M))
kGsto = np.zeros((N,M))
Rin = np.zeros(N)
for m in range(M):
       s = int(Tind[m])
       Gsto[s:N,m] = G0[m]*np.exp(-kS[m]*(t[s:N]-t[s]))
       kGsto[:,m] = kS[m]*Gsto[:,m]
       Rin = Rin + kGsto[:,m]
       GSTO = GSTO + Gsto[:,m]


# GUT ===================================================================
# Gut decay half-life and decay constant
h_gut = 7;       # half-life   [ min]
h_gut = 7/(60*24)
k_gut = -np.log(1/2)/h_gut

Ggut = np.zeros(N)

for c in range(N-1):
     Ggut[c+1] = Ggut[c] + dt*(Rin[c] - k_gut*Ggut[c])

Rgut = k_gut*Ggut


# SOLVE ODEs EULER METHOD  ===============================================
G = np.zeros(N) 
I = np.zeros(N)
B = np.zeros(N)
G[0] = Gss
I[0] = Iss
B[0] = Bss

for c in range(N-1): 
      g1 = k[1]
      g2 = -k[2]*G[c]
      g3 = -k[3]*I[c]*G[c] / (k[4]*G[c]+1)
                   
      G[c+1] = G[c] + dt*(g1 + g2 + g3 + Rgut[c])
     
      g5 = B[c]*k[5]*G[c+1]**2 / (k[6]**2 + G[c+1]**2) 
      g6 = -k[7]*I[c]
      I[c+1] = I[c] + dt*(g5 + g6)
    
      g7  = -k[8]
      g8  = k[9]*G[c+1]
      g9 = -k[10]*G[c+1]**2
      g10 = -k[11]*np.tan(k[12]*B[c])
      B[c+1] = B[c] + dt*( (g7+g8+g9)*B[c] + g10 )
 

# MEAN & LOAD =============================================================
def simpson1d(f,xMin,xMax):
    N = len(f)
    h = (xMax - xMin) / (N - 1)
    
    integral = (h/3) * (f[0] + 2*sum(f[:N-2:2]) \
            + 4*sum(f[1:N-1:2]) + f[N-1])
 
    if N%2 == 0:
        integral = 'N must be an odd number'
        print('integral')
    return integral

Gload = simpson1d(G,0,t2)
Iload = simpson1d(I,0,t2)
Gmean = np.mean(G); Imean = np.mean(I) 

# PEAKS  ===================================================================
peaks, _ = find_peaks(G)
Gpeak = G[peaks]
#tGpeak = t[peaks]*24

Ipeaks, _ = find_peaks(I)
L = len(Ipeaks)
Ipeak = I[Ipeaks]
#tIpeak = t[Ipeaks[L-1]]*24
#Ipeak = I[Ipeaks]
#tIpeak = t[Ipeaks]*24

# hh = 0.5*(Gpeak+Gss)
# cH = np.where(G >= hh)[0][0]
# Gfwhm = dt*(cH[-1] - cH[0] + 1)*24

# hh = 0.5*(Ipeak+Iss)
# cH = np.where(I >= hh)[0]
# Ifwhm = dt*(cH[-1] - cH[0] + 1)*24

# s = int(Tind[0])
# tStart = t[s]*24
# dtG= tGpeak - tStart
# dtI= tIpeak - tStart

# # stomach
# peaks, _ = find_peaks(GSTO)
# GSTOpeak = GSTO[peaks[0]]
# tGSTOpeak = t[peaks[0]]*24
# # gut
# peaks, _ = find_peaks(Ggut)
# Ggutpeak = Ggut[peaks[0]]
# tgutpeak = t[peaks[0]]*24

# time to return to within 5% of basal level
G5 = 1.05*Gss
z = np.where(G>G5)[0][-1]
tG5 = t[z]*24 - T[0]*24
I5 = 1.05*Iss
z = np.where(I>I5)[0][-1]
tI5 = t[z]*24 - T[0]*24

# CONSOLE OUTPUT  =========================================================
print('GLUCOSE')
print('  G_basal = %2.0f' % Gss  )
#print('  G_peak = %2.0f' % Gpeak  )
print('peaks')
print(Gpeak)
print('  G_avg = %2.0f' % Gmean  )
print('  G_load = %2.0f' % Gload  )
# #print('  G_fwhm = %2.2f' % Gfwhm  )
# print('  t_Gpeak = %2.2f' % tGpeak  )
# print('  dt_G = %2.2f' % dtG )
print('  dt_G5 = %2.2f' % tG5 )
print('  ')
print('INSULIN')
print('  I_basal = %2.0f' % Iss  )
#print('  I_peak = %2.0f' % Ipeak  )
print('peaks')
print(Ipeak)
print('  I_avg = %2.0f' % Imean  )
print('  I_load = %2.0f' % Iload  )
# print('  I_fwhm = %2.2f' % Ifwhm  )
# print('  tI_peak = %2.2f' % tIpeak  )
# print('  dt_I = %2.2f' % dtI )
print('  dt_I5 = %2.2f' % tI5 )
# print('  ')
# print('MEAL')
# print('  stomach peak = %2.0f' % GSTOpeak  )
# print('  stomach peak time = %2.2f' % tGSTOpeak  )
# print('  gut = %2.0f' % Ggutpeak  )
# print('  gut peak time = %2.2f' % tgutpeak  )
# print('  dt_I5 = %2.2f' % tI5 )

#%%           
# GRAPHICS ================================================================

#plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.family'] = ['Garamond']
plt.rcParams['font.size'] = 16
font1 = {'color':'black','size':14}

 

fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.23,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  

xLabel = '$t$  [ hours ]'
XLIM = [0,24]
xR = np.arange(0,24.5,2)

# Figure 1:  G vs t ------------------------------------------------------
xP = t*24
yP = G
plt.plot(xP,yP,linewidth=2,color = col)
plt.xlim(XLIM)

plt.xticks(xR)
plt.ylim([50,400])
yR = np.arange(50,405,50) 
plt.yticks(yR)
plt.grid('visible')
plt.ylabel("G   $[ mg.dL^{-1}$ ]", fontdict = font1)
plt.xlabel(xLabel, fontdict = font1)
plt.savefig('a_BIG_53_1.png')

# Figure 2:  I vs t  ------------------------------------------------------
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.23,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  
xP = t*24
yP = I
plt.plot(xP,yP,linewidth=2,color = col)
plt.xlim(XLIM)
plt.xticks(xR)
plt.ylim([10,70])
yR = np.arange(10,71,10) 
plt.yticks(yR)
plt.grid('visible')
plt.ylabel("I   $[ mU.L^{-1}$ ]", fontdict = font1)
plt.xlabel(xLabel, fontdict = font1)
plt.savefig('a_BIG_53_2.png')

# Figure 3:  beta vs t  ------------------------------------------------------
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.23,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  
xP = t*24
yP = B
plt.plot(xP,yP,linewidth=2,color = col)
plt.xlim(XLIM)
plt.xticks(xR)
#plt.ylim([0,100])
#yR = np.arange(50,250,50) 
#plt.yticks(yR)
plt.grid('visible')
plt.ylabel(r"$\beta $  [ mg ]", fontdict = font1)
plt.xlabel(xLabel, fontdict = font1)
plt.savefig('a_BIG_53_3.png')

# Figure 4:  Gsto vs t  ------------------------------------------------------
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.23,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  
xP = t*24
yP = GSTO
plt.plot(xP,yP,linewidth=2,color = col)
plt.xlim(XLIM)
plt.xticks(xR)
plt.ylim([0,600])
yR = np.arange(0,601,100) 
plt.yticks(yR)
plt.grid('visible')
plt.ylabel(r"$G_{sto}$ ", fontdict = font1)
plt.xlabel(xLabel, fontdict = font1)
plt.savefig('a_BIG_53_4.png')

# Figure 5:  Ggut vs t  ------------------------------------------------------
fig = plt.figure(figsize = (4, 3))
fig.subplots_adjust(top=0.94, bottom = 0.20, left = 0.23,\
                    right = 0.96, hspace = 0.60,wspace=0.5)  
xP = t*24
yP = Ggut
plt.plot(xP,yP,linewidth=2,color = col)
plt.xlim(XLIM)
plt.xticks(xR)
plt.ylim([0,100])
yR = np.arange(0,101,10) 
plt.yticks(yR)
plt.grid('visible')
plt.ylabel(r"$G_{gut}$ ", fontdict = font1)
plt.xlabel(xLabel, fontdict = font1)
plt.savefig('a_BIG_53_5.png')