# -*- coding: utf-8 -*-
'''

mnsIZH03.py      FEB 2026

COMPUTATION BIOPHYSICS AND NEUROSCIENCE

IZHIKEVICH MODEL FOR SPIKING NEURAL NETWORK

Ian Cooper
   https://d-arora.github.io/Doing-Physics-With-Matlab/
   
Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/mns/mnsIZH03.pdf

'''

#%% Libraries
import numpy as np
from numpy import pi, sin, cos, linspace, zeros, array, ones, sqrt 
from scipy.integrate import odeint, simpson
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from numpy.linalg import eig
import random
import time

tStart = time.time()
plt.close('all')


#%% INPUTS AND MODEL PARAMETERS
#Ne = 1010; Ni = 240         # number of excitatory & inhibitory neurons
Ne = 1000; Ni = 250
N = Ne+Ni
T = 1001  # imulation time [ms]

# a time scale of the recovery variable u for each neuron and type
# b sensitivity of the recovery variable u to the subthreshold fluctuations of the membrane pot
# c after-spike reset value of the membrane potential for each neuron and type
# d after-spike reset value of the recovery variable u for each neuron and type
# coefficients
a = 0.02*ones((N,1))
b = 0.20*ones((N,1))
c = -65*ones((N,1))
d = 2.0*ones((N,1))

# a = 0.1*ones((N,1))

# Biases random number 0 to 1
np.random.seed(0)
re = np.random.rand(Ne,1)
ri = np.random.rand(Ni,1)

# modified coefficients
a[Ne:N] = 0.02 + 0.08*ri
b[Ne:N] = 0.25 - 0.05*ri
c[0:Ne] = -65 + 15*re**2
d[0:Ne] = 8 - 6*re**2

# Stimulus matrix
S1, S2 = 0.5, -1              # 1
#S1, S2 = 0.60, -1.6         # 2
#S1, S2 = 0.60, -0.6       # 3
#S1, S2 = 0.30, -0.1         # 4
#S1, S2 = 0.10, -0.1         # 5

S = zeros((N,N))
S[:,0:Ne] = S1 * np.random.rand(Ne+Ni, Ne)
S[:,Ne:N] = S2*np.random.rand(Ne+Ni, Ni)

# Initialize arrays
v = -65 * np.ones((Ne+Ni, 1))
u = b * v
I = zeros((N,1))
firings = zeros((N,2))      # spiking neurons
I_array = np.zeros((N, T))
v_array = np.zeros((N, T))
u_array = np.zeros((N, T))

fired = np.where(v >= 30)[0]   # indices of fired neurons
spikes = zeros(T)              # count number of fired neurons each time step

# Time evolution of network:  1000 ms
for t in range(0, T):
    # Step 1: input current:
    #   calculate input current for each neuron with noise contribution I_ext(t)
    I[0:Ne] = 5 * np.random.randn(Ne, 1)
    I[Ne:N] = 2 * np.random.randn(Ni, 1)
    # Sum synaptic contributions for fired neurons in previous time step
    if t > 0:  
       I += np.sum(S[:, fired], axis=1).reshape(-1, 1)
      
    # Step 2: update the membrane pot and recovery variable with Euler's method
    v += 0.5 * (0.04 * v**2 + 5 * v + 140 - u + I)
    v += 0.5 * (0.04 * v**2 + 5 * v + 140 - u + I)
    u +=  a * (b * v - u)
    
    # Step 3: check for spikes and update the membrane pot and recovery variable:
    fired = np.where(v >= 30)[0] # check if the membrane potential exceeds 30 mV
    spikes[t] = len(fired) #np.sum(v >= 30)
    if fired.size > 0:
        firings = np.vstack((firings, np.hstack((t * np.ones((fired.size, 1)), fired.reshape(-1, 1)))))
    # equalize all spikes at 30 mV by resetting v first to +30 mV and then to c:
    v[fired] = c[fired] # reset v for fired neurons
    u[fired] = u[fired] + d[fired] # increment u for fired neurons
    
    v = np.clip(v, -100, 30) 
    
    # Step 4: Record data for each time step
    I_array[:, t] = I.flatten()
    v_array[:, t] = v.flatten()
    u_array[:, t] = u.flatten()


#%% Fourier Transformand psd
Nf = T
f = linspace(1,60,Nf)
tF = linspace(0,999/1001,Nf)
H = zeros((Nf,1))
h = sin(2*pi*2*tF)
for cc in range(Nf):
     g = spikes * np.exp(1j*2*pi*f[cc]*tF)
     H[cc] = simpson(np.real(g),tF) + 1j*simpson(np.imag(g),tF)

psd = np.conj(H)*H
    
#%%   FIG 1: fired neurons as a function of time 
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)

ax.scatter(firings[:, 0], firings[:, 1], s=3, c='b')
ax.scatter(firings[firings[:,1]>=Ne, 0], firings[firings[:,1]>=Ne,1], s=1, c='r')
ax.scatter(firings[firings[:,1]<Ne, 0], firings[firings[:,1]<Ne,1], s=1, c='b')
ax.set_xlabel('t  [ ms ]')
ax.set_ylabel('neuron number')
ax.axhline(y = Ne, color='k', linestyle='-', linewidth=1)
ax.text(750, 900,'EXCITATORY',bbox=dict(facecolor='white', alpha=1))
ax.text(750, 1100,'INHIBITORY',bbox=dict(facecolor='white', alpha=1))
fig1.tight_layout()
fig1.savefig('a1.png')


#%%   FIG 2: Percentage spiking neurons   
plt.rcParams["figure.figsize"] = (6,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)

ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('% spiking neurons')
xP = range(T); yP = spikes/10
ax.plot(xP,yP,'b',lw = 1)
fig2.tight_layout()
fig2.savefig('a2.png')


#%%   FIG 3: time evolution v   
plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=2, ncols=1)
ax[0].set_title('membrane potential v',fontsize = 11)
xP = range(T); yP = v_array[100,:]
ax[0].plot(xP,yP,'b', lw = 1, label ='excitatory')
ax[0].legend(fontsize = 10)
ax[1].set_xlabel('t  [ ms ]'); ax[1].set_ylabel('v  [ mV ]')
yP = v_array[1100,:]
ax[1].plot(xP,yP,'r', lw = 1, label = 'inhibitory')
ax[1].legend(fontsize = 10)
fig3.tight_layout()
fig3.savefig('a3.png')


#%%   FIG 3: time evolution  u  
plt.rcParams["figure.figsize"] = (6,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_title('recovery variable u',fontsize = 11)
ax.set_xlabel('t  [ ms ]'); ax.set_ylabel('u')
xP = range(T); yP = u_array[100,:]
ax.plot(xP,yP,'b', lw = 1, label ='excitatory')
yP = u_array[1100,:]
ax.plot(xP,yP,'r', lw = 1, label = 'inhibitory')
ax.legend(fontsize = 10)
fig4.tight_layout()
fig4.savefig('a4.png')


#%% FIG 5: exteral current stimulus   
plt.rcParams["figure.figsize"] = (6,3)
fig5, ax = plt.subplots(nrows=1, ncols=1)
plt.imshow(I_array, aspect='auto', cmap='seismic', origin='lower') 
ax.set_xlabel('t  [ms]')
ax.set_ylabel('neuron number')
plt.colorbar(label='current I (pA)')
ax.set_title('Current stimulus',fontsize = 12)
fig5.tight_layout()
fig5.savefig('a5.png')

#%% psd
plt.rcParams["figure.figsize"] = (6,3)
fig6, ax = plt.subplots(nrows=1, ncols=1)
ax.set_title('Foureir Transform',fontsize = 11)
ax.set_xlabel('f  [Hz]'); ax.set_ylabel('psd')
ax.grid()
xP = f; yP = psd
ax.plot(xP,yP,'b', lw = 1)
fig6.tight_layout()
fig6.savefig('a6.png')

         
#%%
tExe = time.time() - tStart
print('  ')
print('Execution time')
print(tExe)


