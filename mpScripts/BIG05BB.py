# BIG05.py
# -*- coding: utf-8 -*-
# Ian Cooper
# 12 Nov 2023
# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Beta-cell mass B, insulin I, glucose G dynamics
# Time evolution of beta-cell mass, insulin, gluocose
# Solution of 3 coupled ODEs using the odeint function

# Response of glucose regulatory system after the consumption of food

# LIBRARIES  ==============================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.signal import find_peaks

plt.show(block = True)

# MEALS (STOMACH) ======================================================
# Inputs for 8 meals: Quantity G0 / Half life h > 0 / Meal time T
# >>> Meal quantity: glucose
G0 = np.array([300,00,0,0,0,0,0,0])  

# >>> Half-life [min] h > 0   
h  = np.array([30,20,1,1,1,1,1,1,])  

# >>> Time of meal  [hours]    
T  = np.array([2,10,0,0,0,0,0,0,])   

M = 8    # max number of meals

# SETUP: TIME SCALE / MEALS ================================================
# TIME SCALE [minutes]
nT =  3600*24+1              # Grid points: number of sec in a day + 1 
tMax = 24*60                 # Number of minutes in a day
t = np.linspace(0,tMax,nT)   # time span [min]
dt = t[2] - t[1]             # time increnebt [ min ]
f = 1/(24*60)                # 1 / (minutes in days)

k = np.log(2)/h              # Stomach: decay constant [1/min]
R  = np.zeros(M)             # Index for start time
R = T*3600                   # convert hours to min

# 8 MEALS / m meal number
GSTO = np.zeros(nT)
Gsto = np.zeros((nT,M))
kGsto = np.zeros((nT,M))
Rin = np.zeros(nT)
for m in range(M):
       s = int(R[m])       # Index for start time
       Gsto[s:nT,m] = G0[m]*np.exp(-k[m]*(t[s:nT]-t[s]))
       kGsto[:,m] = k[m]*Gsto[:,m]
       Rin = Rin + kGsto[:,m]
       GSTO= GSTO + Gsto[:,m]
                          
# Gut decay half-life and decay constant
h_gut = 7;       # half-life   [ day]
k_gut = np.log(2)/h_gut

Ggut = np.zeros(nT)

for c in range(nT-1):
     Ggut[c+1] = Ggut[c] + dt*(Rin[c] - k_gut*Ggut[c])

Rgut = k_gut*Ggut


# SOLVE odes  ============================================================

# FUNCTIONS  ==============================================================
#beta-cell mass (z) insulin (y)   glucose (x)  
def lorenz(t, state, k,):    
    x, y, z = state
    global c, Rgut
    c = c + 1
    
    dx = k[1] - k[2]*x - k[3]*x*y/(k[4]*x+1) + Rgut[c-1] 
    dy = k[5]*z*x**2/(k[6]+x**2) - k[7]*y
    dz = (-k[8] + k[9]*x - k[10]*x**2 )*z - k[11]*np.tan(k[12]*z)
    
    print(c)
    return [dx, dy, dz]


# CONSTANTS  ===============================================================
k = np.zeros(13)
k[1]  = 864*f   # 864
k[2]  = 1.44*f     # 1.44
k[3]  = 1*f # 1
k[4]  = 0.01     # 0.01
k[5]  = 40*f # 40
k[6]  = 20000    # 20000
k[7]  = 432*f      # 432
k[8]  = 0.06*f    # 0.06
k[9]  = 10e-4*f  # 10e-4
k[10] = 2.4e-6*f   # 2.4e-6
k[11] = 3.49*f    # 3.49
k[12] = np.pi/3050    # pi/3050


#%%
# FIXED POINTS: Steady values Gss, Iss, Bss 
# Must select starting values carefully to obtain correct solution

def myFunction(u):
   x = u[0]
   y = u[1]
   z = u[2]
   F = np.empty((3))
   F[0] = k[1] - k[2]*x - k[3]*y*x/(k[4]*x+1)
   F[1] = z*k[5]*x**2 / (k[6]+x**2)-k[7]*y
   F[2] = (-k[8] + k[9]*x - k[10]*x**2)*z - k[11]*np.tan(k[12]*z)
   return F

# >>> Statring values Gss, Iss, Bss
zGuess = np.array([80,17,755])     

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


#%% SOLVE ODEs EULER METHOD
G = np.zeros(nT) 
I = np.zeros(nT)
B = np.zeros(nT)
G[0] = Gss
I[0] = Iss
B[0] = Bss

for c in range(nT-1): 
      g1 = k[1]
      g2 = -k[2]*G[c]
      g3 = -k[3]*I[c]*G[c] / (k[4]*G[c]+1)
      g4 = k_gut*Ggut[c]
      G[c+1] = G[c] + dt*(g1 + g2 + g3 + g4)
     
      g5 = B[c]*k[5]*G[c+1]**2 / (k[6] + G[c+1]**2) 
      g6 = -k[7]*I[c]
      I[c+1] = I[c] + dt*(g5 + g6)
    
      g7  = -k[8]
      g8  = k[9]*G[c+1]
      g9 = -k[10]*G[c+1]**2
      g10 = -k[11]*np.tan(k[12]*B[c])
      B[c+1] = B[c] + dt*( (g7+g8+g9)*B[c] + g10 )
      

#%% PEAK VALUES: Glucose and Insulin 

indG = find_peaks(G)[0]
indI = find_peaks(I)[0]
Gpeaks = G[indG]                   # peaks in glucose
tGpeaks = t[indG]/60               # times for glucose peaks  [hours]
Ipeaks = I[indI]                   # peaks in insulin [hours]
tIpeaks = t[indI]/60               # times for insulin peaks [hours]
print('Peaks')
print('Gpeaks = %5.2f' % Gpeaks[0], 'tGpeaks = %5.2f' % tGpeaks[0]  )
print('Ipeaks = %5.2f' % Ipeaks[0], 'tIpeaks = %5.2f' % tIpeaks[0]  )


#%% GLUOCSE AND INSULIN LOADS / MEANS   
def simpson1d(f,xMin,xMax):
    N = len(f)
    h = (xMax - xMin) / (N - 1)
    
    integral = (h/3) * (f[0] + 2*sum(f[:N-2:2]) \
            + 4*sum(f[1:N-1:2]) + f[N-1])
 
    if N%2 == 0:
        integral = 'N must be an odd number'
        print('integral')
    return integral

fG = G; fI = I
Gload = simpson1d(fG,0,24)
Iload = simpson1d(fI,0,24)
Gmean = np.mean(G); Imean = np.mean(I)

print('  ')
print('Gload = %5.2f ' % Gload, 'Iload = %5.2f ' % Iload ) 
print('Gmean = %5.2f ' % Gmean, 'Imean = %5.2f ' % Imean ) 

#%%  GRAPHICS  

def plotFn(xP,yP,limX,tX,limY,tY,txtx,txtY):
    font1 = {'family':'Cambria','color':'black','size':12}
    plt.rcParams['font.family'] = ['Tahoma']
    plt.rcParams['font.size'] = 12

    fig = plt.figure(figsize=(5,3))
    fig.subplots_adjust(top=0.90, bottom = 0.20, left = 0.20,\
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
    
#%%   STOMACH: MEAL INPUT
y1 = 0; y2 = 510; dy = 100
x1 = 0; x2 = 25; dx = 2
limX = [x1, x2]
tX = range(x1,x2,dx)
limY= [y1, y2]
tY = range(y1,y2,dy)
txtX = '$t$    [hours]'
txtY = '$G_{sto}$    $[ mg.dL^{-1} ]$ '
xP = t/60
yP = GSTO
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)

#%%   GUT: MEAL INPUT
y1 = 0; y2 = 102; dy = 10
x1 = 0; x2 = 25; dx = 2
limX = [x1, x2]
tX = range(x1,x2,dx)
limY= [y1, y2]
tY = range(y1,y2,dy)
txtX = '$t$    [hours]'
txtY = '$G_{gut}$    $[ mg.dL^{-1} ]$ '
xP = t/60
yP = Ggut
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)

#%%   PLASMA GLUCOSE
y1 = 0; y2 = 550; dy = 100
x1 = 0; x2 = 25; dx = 2
limX = [x1, x2]
tX = range(x1,x2,dx)
limY= [y1, y2]
tY = range(y1,y2,dy)
txtX = '$t$    [hours]'
txtY = '$G$    $[ mg.dL^{-1} ]$ '
xP = t/60
yP = G
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)

#%%   PLASMA INSULIN
y1 = 0; y2 = 105; dy = 20
x1 = 0; x2 = 25; dx = 2
limX = [x1, x2]
tX = range(x1,x2,dx)
limY= [y1, y2]
tY = range(y1,y2,dy)
txtX = '$t$    [hours]'
txtY = '$I$    $[ mU.L^{-1} ]$ '
xP = t/60
yP = I
plotFn(xP,yP,limX,tX,limY,tY,txtX,txtY)

# # GRAPHICS  ===============================================================
# #font1 = {'family':'Tahoma','color':'black','size':12}
# font1 = {'family':'Cambria','color':'black','size':12}
# #font2 = {'family':'serif','color':'blue','weight':'normal','size':14}
# #plt.rcParams['font.family'] = ['Cambria']
# plt.rcParams['font.family'] = ['Tahoma']
# plt.rcParams['font.size'] = 12

# fig2 = plt.figure(figsize=(5,5))
# fig2.subplots_adjust(top=0.90, bottom = 0.20, left = 0.20,\
#                     right = 0.94, hspace = 0.45,wspace=0.5)

# ax = fig2.add_subplot(2, 1, 1)    
# xP = t/60
# yP = G

# plt.plot(xP,yP,linewidth=2,color='b')


# #plt.xlim([0,40])

# #plt.ylim([0,300])
# plt.grid('visible')
   
# plt.xticks(np.arange(0,25,2))
# #plt.legend()
# plt.xlabel('t [ days ]',fontdict = font1)
# #plt.ylabel('$G_{meal}$ ', fontdict = font1)
# plt.ylabel('G ', fontdict = font1)
# plt.title('  ')
# plt.show()



# ax = fig2.add_subplot(2, 1, 2)    
# xP = t/60
# yP = I

# plt.plot(xP,yP,linewidth=2,color='b')

# #plt.xlim([0,40])

# #plt.ylim([0,300])
# plt.grid('visible')
   
# plt.xticks(np.arange(0,25,2))
# #plt.legend()
# plt.xlabel('t [ days ]',fontdict = font1)
# #plt.ylabel('$G_{meal}$ ', fontdict = font1)
# plt.ylabel('I ', fontdict = font1)
# plt.title('  ')
# plt.show()