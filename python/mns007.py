# -*- coding: utf-8 -*-
'''
mns007.py      Jan 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE
   LEAK / FAST Na+ Channels
GEOMETRICAL ANAYLSIS OF [1D] DYNAMICAL SYSTEMS

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mns007.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt 
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.linalg import eig
import time

tStart = time.time()
plt.close('all')


#%% FUNCTIONS  Solve ODE for x,y    
def lorenz(t, state): 
    x = state
    dx = -( GNa*(x - ENa) + GL*(x - EL) + I0 ) / CM
    return dx  

    
#%% Time span
s = 1e3    # Change of units  
N = 9999
tS = 1e-3
t = linspace(0,tS,N)
h = t[2] - t[1]

# Initial membrane voltage [V]  / external current [A]   I0 == Iext
V0 = -100e-3
I0 = -0.6e-3

# Model parameters
CM = 10e-6      # Membrane capacitance [F]
GL = 19e-3      # leak conductance [S]
GNa = 74e-3    # Na+ conductance  [S]
EL = -67e-3     # Leak reversal potential / Nernst potential  [V]
ENa = 60e-3     # Na+ reversal potential

#%%  SETUP
Vss = (-I0 + GNa*ENa + GL*EL) / (GNa + GL)

IL_ss = GL*(Vss - EL)
INa_ss = GNa*((Vss - ENa))
IM_ss = IL_ss + INa_ss
IC_ss = IM_ss + I0

Vdot = -( GNa*(Vss - ENa) + GL*(Vss - EL) + I0 ) / CM

q = Vss*s; print('Vss = %0.1f mV' %q)
q = IC_ss*s; print('IC_ss = %0.1f mA' %q)
q = INa_ss*s; print('INa_ss = %0.1f mA' %q)
q = IL_ss*s; print('IL_ss = %0.1f mA' %q)
q = IM_ss*s; print('IM_ss = %0.1f mA' %q)
q = I0*s; print('I0 = %0.1f mA' %q)

#%% Solve ODE
sol = odeint(lorenz, V0, t, tfirst=True)
VM = sol[:,0]     

IL =  GL*(VM - EL)
INa = GNa*(VM - ENa)
Iext = ones(N)*I0
IC = -(INa + IL + Iext)
IM = INa + IL

#%%   FIG 1: t vs VM   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms ]')
ax.set_ylabel('V$_M$ [mV]')
q = I0*1e3;ax.set_title('I$_{ext}$ = %0.2f mA' %q)
q = ENa*s; ax.text(0.7,70,'E$_{Na}$ = %0.1f mV' %q)
q = EL*s; ax.text(0.7,-50,'E$_L$ = %0.1f mV' %q)
ax.grid()
ax.set_ylim([-100,100])
xP = t*s; yP = VM*s; ax.plot(xP, yP,'b',lw = 2)
xP = [0,tS*s]; yP = [EL*s,EL*s]; ax.plot(xP,yP,'r',lw = 1)
xP = [0,tS*s]; yP = [ENa*s,ENa*s]; ax.plot(xP,yP,'r',lw = 1)
fig1.tight_layout()
fig1.savefig('a1.png')

#%%   FIG 2: t vs I
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5.5,3.5)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('t  [ ms ]')
ax.set_ylabel('I  [ mA ]')
q1 = I0*s; q2 = V0*s
ax.set_title('I$_{ext}$ = %0.2f mA' %q1 + '   V$_M$(0) = %0.f mV' %q2,
             fontsize = 10)
ax.grid()

xP = t*s; yP = Iext*s; ax.plot(xP, yP,'k',lw = 2,label = 'Iext')
xP = t*s; yP = IL*s; ax.plot(xP, yP,'r',lw = 2,label = 'IL')
xP = t*s; yP = INa*s; ax.plot(xP, yP,'b',lw = 2,label = 'INa')
xP = t*s; yP = IC*s; ax.plot(xP, yP,'g',lw = 2,label = 'IC')
xP = t*s; yP = IM*s; ax.plot(xP, yP,'m',lw = 2,label = 'IM')

ax.legend( loc='upper center', 
    bbox_to_anchor=(0.5, 1.35),
    fontsize = 10,
    ncol=5, )
fig2.tight_layout()
fig2.savefig('a2.png')

#%%   FIG 3:  Phase portrait
V1 = -50e-3; V2 = 100e-3
V = linspace(V1,V2,N)
vDot = -( GNa*(V - ENa) + GL*(V - EL) + I0 ) / CM

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,2.5)
fig3, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel(r'v$_{dot}$  [ mV.ms$^{-1}$ ]')
ax.set_xlabel('V$_M$ [ mV ]')
q = I0*1e3;ax.set_title('I$_{ext}$ = %0.2f mA' %q)
ax.grid()
#ax.set_ylim([-130,70])
ax.plot(V*s,vDot,'b',lw = 2)
xP = V[vDot<0][0]*1e3; yP = 0 
ax.plot(xP,yP,'ro',ms = 8)
fig3.tight_layout()
fig3.savefig('a3.png')

#%%   FIG 4: I-V characteristic I vs V   
V1 = -100e-3; V2 = 200e-3
V = linspace(V1,V2,N)
IINa = GNa*(V - ENa)
IIL =  GL*(V - EL)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (5,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)

ax.set_ylabel(r'I  [ mA ]')
ax.set_xlabel('V$_M$ [ mV ]')
ax.grid()
ax.set_ylim([-10,10]); ax.set_xlim([-100,200])
ax.plot(V*s,IIL*s,'r',lw = 2,label = 'I$_L$')
ax.plot(V*s,IINa*s,'b',lw = 2,label = 'I$_{Na}$')
ax.plot([Vss*s,Vss*s],[-10,10],'g',lw = 1,label = 'V$_{SS}$')
ax.plot([EL*s,EL*s],[-10,10],'r',lw = 1,  label = 'V$_L$')
ax.plot([ENa*s,ENa*s],[-10,10],'b',lw = 1,label = 'V$_{Na}$')
ax.legend(frameon = False)
fig4.tight_layout()
fig4.savefig('a4.png')


#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


