# -*- coding: utf-8 -*-
# BIG02.py
# Ian Cooper
# 5 Nov 2023
# Website: https://d-arora.github.io/Doing-Physics-With-Matlab/

# Beta-cell mass B, insulin I, glucose G dynamics
# Time evolution of beta-cell mass, insulin, gluocose
# Phase plot  G vs I at a constant beta-cell mass
# Solution of 3 coupled ODEs using the odeint function


# LIBRARIES  ==============================================================
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
 

# FUNCTIONS  ==============================================================
#beta-cell mass (z) insulin (y)   glucose (x)  
def lorenz(t, state, k):    
    x, y, z = state
    
    dx = k[1] - k[2]*x - k[3]*x*y/(k[4]*x+1) 
    dy = k[5]*z*x**2/(k[6]**2+x**2) - k[7]*y
    dz = (-k[8] + k[9]*x - k[10]*x**2 )*z - k[11]*np.tan(k[12]*z)

# # k10    
#     u = k[10]      
#     if t > 50:
#         u = 3.4e-6
#     dx = k[1] - k[2]*x - k[3]*x*y/(k[4]*x+1) 
#     dy = k[5]*z*x**2/(k[6]+x**2) - k[7]*y
#     dz = (-k[8] + k[9]*x - u*x**2)*z

# k3    insulin reistance    
    # w = k[3]      
    # if t > 50:
    #     w = 0.3
    # dx = k[1] - k[2]*x - w*x*y/(k[4]*x+1) 
    # dy = k[5]*z*x**2/(k[6]+x**2) - k[7]*y
    # dz = (-k[8] + k[9]*x - k[10]*x**2 )*z - k[11]*np.tan(k[12]*z)
    
    
# k5    insulin secretion     
    # w = k[5]      
    # if t > 50:
    #     w = 40
    # dx = k[1] - k[2]*x - k[3]*x*y/(k[4]*x+1) 
    # dy = w*z*x**2/(k[6]+x**2) - k[7]*y
    # dz = (-k[8] + k[9]*x - k[10]*x**2 )*z - k[11]*np.tan(k[12]*z)    
    

# k3 and k5    insulum  resistance / insulin secretion     
    # w = k[3]  
    # u = k[5]    
    # if t > 20:
    #     w = 0.3
    # if t > 120:
    #     u = 40
    
    # dx = k[1] - k[2]*x - w*x*y/(k[4]*x+1) 
    # dy = u*z*x**2/(k[6]+x**2) - k[7]*y
    # dz = (-k[8] + k[9]*x - k[10]*x**2 )*z - k[11]*np.tan(k[12]*z)    
           
    
    return [dx, dy, dz]

 #%%   
# CONSTANTS  ===============================================================
k = np.zeros(13)
k[1]  = 864      # 864
k[2]  = 1.44     # 1.44
k[3]  = 1        # 1
k[4]  = 0.01     # 0.01
k[5]  = 40       # 40
k[6]  = 150      # 150
k[7]  = 432      # 432
k[8]  = 0.06     # 0.06
k[9]  = 10e-4    # 10e-4
k[10] = 2.4e-6   # 2.4e-6
k[11] = 3.49     # 3.49
k[12] = np.pi/3050    # pi/3050


# SOLVE ODEs  ============================================================ 
#y0 = [80, 17, 755]
y0 = [80, 17, 755]  # Inital Conditions: G0  I0   B0
N = 999            # Time interval for t
t1 = 0
t2 = 500
tSpan = np.linspace(t1,t2,N)
t = tSpan

sol = odeint(lorenz, y0, tSpan, args = (k,), tfirst=True)

G = sol[:,0]
I = sol[:,1]
B = sol[:,2]


print('  ')
print('Physiological fixed-point')
print('Bss = %5.2f ' % B[N-1],  'Iss = %5.2f' % I[N-1], 'Gss = %5.2f' % G[N-1]  )
print('  ')


# PHASE PLOT nullclines and critical point (Ic, Gc)======================
Gn = np.linspace(10,600,599);
IG = (k[1]-k[2]*Gn)*(k[4]*Gn+1)/(0.5*k[3]*Gn);    # G nullcine   k[3] ??
II = B[N-1]*k[5]*Gn**2/((k[6] + Gn**2)*k[7])  # I nullcline
D =IG-II
z = np.where(D<0)[0][0]
Gc = Gn[z]       # critical point
Ic = II[z]        
  
#%%           
# GRAPHICS ================================================================
#plt.close('all')
font1 = {'family':'Tahoma','color':'black','size':12}
plt.rcParams['font.family'] = ['Tahoma']
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize = (4, 4.6))
fig.subplots_adjust(top=0.94, bottom = 0.16, left = 0.23,\
                    right = 0.96, hspace = 0.60,wspace=0.5)

ax = fig.add_subplot(3, 1, 1)
plt.plot(t,G,linewidth=2,color='b')
plt.ylim([0,250])
yR = np.arange(0,250,100) 
plt.yticks(yR)
plt.grid('visible')
plt.ylabel("G   $[ mg.dL^{-1}$ ]", fontdict = font1)
#plt.xlabel('t [days]', fontdict = font1)

ax = fig.add_subplot(3, 1, 2)
plt.plot(t,I,linewidth=2,color='b')
plt.ylim([0,50])
yR = np.arange(0,60,20) 
plt.yticks(yR)
plt.grid('visible')
#plt.xlabel('t   $[days]$', fontdict = font1)
plt.ylabel("I   $[ mU.L^{-1}$ ]", fontdict = font1)

ax = fig.add_subplot(3, 1, 3)
plt.plot(t,B,linewidth=2,color='b')
plt.ylim([0,1500])
yR = np.arange(0,2000,500) 
plt.yticks(yR)
plt.xlabel('t   $[days]$', fontdict = font1)
plt.ylabel("$\\beta$  [ mg ]", fontdict = font1)
plt.grid('visible')
plt.show()


#  GRAPHICS: Phase plot,Nullclines, critical point (Ic, Gc) =============
fig2 = plt.figure(figsize=(4,4))
fig2.subplots_adjust(top=0.88, bottom = 0.20, left = 0.26,\
                    right = 0.94, hspace = 0.45,wspace=0.5)
plt.plot(IG,Gn,linewidth=2,color='m',label='G')
plt.plot(II,Gn,linewidth=2,color='r',label='I')
plt.plot(I,G,linewidth=2,color='b',label='GI')
plt.xlim([0,50])
plt.ylim([0,300])
#xR = [0, 10, 20 , 30, 40, 50]
xR = np.arange(0,60,10)
plt.xticks(xR,fontsize = '12')
yR = np.arange(0,350,50)
plt.yticks(yR,fontsize = '12')
plt.grid('visible')
plt.legend(fontsize = '10', ncol=3)
plt.xlabel('I $[ mU.L^{-1}$ ]',fontdict = font1,fontsize = '12')
plt.ylabel("G   $[ mg.dL^{-1}$ ]", fontdict = font1)
plt.title(f'$\\beta$ = {B[N-1]:2.0f}   $G_c$ = {Gc:2.0f} \
          $I_c$ ={Ic:2.0f}',fontsize=12)
plt.show()

plt.title('  ')
