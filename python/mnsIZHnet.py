# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 20:09:58 2026

@author: Owner
"""

"""
A Spiking Neural Network (SNN) model based on the Izhikevich model (2003).

author: Fabrizio Musacchio
date: Apr 20, 2024
"""



# %% IMPORTS
import numpy as np
import matplotlib.pyplot as plt
# set global font size for plots:
plt.rcParams.update({'font.size': 12})
# create a folder "figures" to save the plots (if it does not exist):
import os
if not os.path.exists('figures'):
    os.makedirs('figures')
import time
    
#%%
tStart = time.time()
plt.close('all')    
# %% MAIN
# for reproducibility:
np.random.seed(0) #100

# simulation time:
T = 1000  # ms

# constants:
Ne = 800  # Number of excitatory neurons
Ni = 202  # Number of inhibitory neurons

# initialize parameters; pre-define parameters for different neuron types:
re = np.random.rand(Ne, 1) # excitatory neurons, "r" stands for random
ri = np.random.rand(Ni, 1) # inhibitory neurons
p_RS  = [0.02, 0.2, -65, 8, "regular spiking (RS)"] # regular spiking settings for excitatory neurons (RS)
p_IB  = [0.02, 0.2, -55, 4, "intrinsically bursting (IB)"] # intrinsically bursting (IB)
p_CH  = [0.02, 0.2, -51, 2, "chattering (CH)"] # chattering (CH)
p_FS  = [0.1, 0.2, -65, 2, "fast spiking (FS)"] # fast spiking (FS)
p_TC  = [0.02, 0.25, -65, 0.05, "thalamic-cortical (TC)"] # thalamic-cortical (TC) (doesn't work well)
p_LTS = [0.02, 0.25, -65, 2, "low-threshold spiking (LTS)"] # low-threshold spiking (LTS)
p_RZ  = [0.1, 0.26, -65, 2, "resonator (RZ)"] # resonator (RZ)
a_e, b_e, c_e, d_e, name_e = p_RS
a_i, b_i, c_i, d_i, name_i = p_LTS
a = np.vstack((a_e * np.ones((Ne, 1)), a_i + 0.08 * ri))
b = np.vstack((b_e * np.ones((Ne, 1)), b_i - 0.05 * ri))
c = np.vstack((c_e + 15 * re**2,       c_i * np.ones((Ni, 1))))
d = np.vstack((d_e-6 * re**2,          d_i * np.ones((Ni, 1))))
# define various synaptic weights for different scenarios (uncomment one of the following):
S = np.hstack((0.5 * np.random.rand(Ne+Ni, Ne), -1*np.random.rand(Ne+Ni, Ni))) # default values from Izhikevich
#S = np.hstack((0.60 * np.random.rand(Ne+Ni, Ne), -1.6*np.random.rand(Ne+Ni, Ni))) # S1
#S = np.hstack((0.60 * np.random.rand(Ne+Ni, Ne), -0.6*np.random.rand(Ne+Ni, Ni))) # S2
#S = np.hstack((0.30 * np.random.rand(Ne+Ni, Ne), -0.1*np.random.rand(Ne+Ni, Ni))) # S3
#S = np.hstack((0.10 * np.random.rand(Ne+Ni, Ne), -0.1*np.random.rand(Ne+Ni, Ni))) # S4
#S = np.hstack((0.30 * np.random.rand(Ne+Ni, Ne), -0.2*np.random.rand(Ne+Ni, Ni))) # S5
"""default values for chattering neurons (regular spiking):
a is the time scale of the recovery variable u for each neuron and type
b is the sensitivity of the recovery variable u to the subthreshold fluctuations of the membrane potential
c is the after-spike reset value of the membrane potential for each neuron and type
d is the after-spike reset value of the recovery variable u for each neuron and type

parameters from the original Izhikevich model for chattering:
a = np.vstack((0.02 * np.ones((Ne, 1)), 0.02 + 0.08 * ri))
b = np.vstack((0.2 * np.ones((Ne, 1)),  0.25 - 0.05 * ri))
c = np.vstack((-65 + 15 * re**2,       -65 * np.ones((Ni, 1))))
d = np.vstack((8 - 6 * re**2,            2 * np.ones((Ni, 1))))
"""


# initial values of v and u:
v = -65 * np.ones((Ne+Ni, 1))
u = b * v
firings = np.array([]).reshape(0, 2)  # Spike timings



# initialize variables for recording data:
I_array = np.zeros((Ne+Ni, T))
v_array = np.zeros((Ne+Ni, T))
u_array = np.zeros((Ne+Ni, T))

fired = 0

# simulation of 1000 ms:
for t in range(0, T):
    # step 1: input current calculation:
    # i.e., calculate input current for each neuron with noise contribution (this is our "I_external(t)").
    I = np.vstack((5 * np.random.randn(Ne, 1), 2 * np.random.randn(Ni, 1)))
    # summing synaptic contributions if there are any fired neurons in previous time step:
    if t > 0:  
        I += np.sum(S[:, fired], axis=1).reshape(-1, 1)
       
    # step 2: update the membrane potential and recovery variable (neuron dynamics) with Euler's method:
    v += 0.5 * (0.04 * v**2 + 5 * v + 140 - u + I)
    v += 0.5 * (0.04 * v**2 + 5 * v + 140 - u + I)
    u +=  a * (b * v - u)
    
    # step 3: check for spikes and update the membrane potential and recovery variable:
    fired = np.where(v >= 30)[0] # check if the membrane potential exceeds 30 mV
    if fired.size > 0:
        firings = np.vstack((firings, np.hstack((t * np.ones((fired.size, 1)), fired.reshape(-1, 1)))))
        # equalize all spikes at 30 mV by resetting v first to +30 mV and then to c:
    v[fired] = c[fired] # reset v for fired neurons
    u[fired] = u[fired] + d[fired] # increment u for fired neurons
    
    """ # clip v to +30 mV if it exceeds this value:
    v = np.clip(v, -100, 30) """
    
    # step 4: record data:
    I_array[:, t] = I.flatten()
    v_array[:, t] = v.flatten()
    u_array[:, t] = u.flatten()

# plots:
a_str = f"{a[0, 0].round(2)}_{a[-1:, 0][0].round(2)}"
b_str = f"{b[0, 0].round(2)}_{b[-1:, 0][0].round(2)}"
c_str = f"{c[0, 0].round(2)}_{c[-1:, 0][0].round(2)}"
d_str = f"{d[0, 0].round(2)}_{d[-1:, 0][0].round(2)}"

# plotting the spike timings:
plt.figure(figsize=(5, 5))
plt.scatter(firings[:, 0], firings[:, 1], s=1, c='k')
excitatory = firings[:, 1] < Ne
inhibitory = firings[:, 1] >= Ne
""" plt.scatter(firings[excitatory, 0], firings[excitatory, 1], s=1, c='blue', label='Excitatory')
plt.scatter(firings[inhibitory, 0], firings[inhibitory, 1], s=1, c='red', label='Inhibitory') """
# plot a horizontal line to separate excitatory and inhibitory neurons:
plt.axhline(y=Ne, color='k', linestyle='-', linewidth=1)
# indicate excitatory and inhibitory neurons with text annotations with background color white:
plt.text(0.8, 0.76, 'excitatory', color='k', fontsize=12, ha='left', va='center', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=1)) 
plt.text(0.8, 0.84, 'inhibitory', color='k', fontsize=12, ha='left', va='center', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=1))

plt.xlabel('Time (ms)')
plt.ylabel('Neuron index')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlim([0, T])
plt.ylim([0, Ne+Ni])
plt.yticks(np.arange(0, Ne+Ni+1, 200))
#plt.legend(markerscale=10, loc='upper right')
plt.tight_layout()
#plt.savefig(f'figures/izhikevich_SNN_firings_Ne_{name_e}_Ni_{name_i}.png', dpi=120)
plt.show()

# %% MORE PLOTS
# plotting the applied current at each neuron and time point:
plt.figure(figsize=(5, 5))
plt.imshow(I_array, aspect='auto', cmap='seismic', origin='lower')
# make the colorbar limits symmetric by choosing the maximum absolute value:
plt.clim(-np.max(np.abs(I_array)), np.max(np.abs(I_array)))
plt.xlabel('Time (ms)')
plt.ylabel('Neuron index')
plt.colorbar(label='Current (pA)')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlim([0, T])
plt.ylim([0, Ne+Ni])
plt.yticks(np.arange(0, Ne+Ni+1, 200))
plt.tight_layout()
plt.savefig(f'figures/izhikevich_SNN_currents_Ne_{name_e}_Ni_{name_i}.png', dpi=120)
plt.show()

# plot the membrane potential v(t) and recovery variable u(t) for one excitatory neuron vs. time:
neuron_e_id = 1
neuron_i_id = 811
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(v_array[neuron_e_id, :], label=f'excitatory neuron (id {neuron_e_id})', color='blue')
plt.title(f'Membrane potential v(t) for typical example neurons')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlim([0, T])
plt.ylabel('$v(t)$ [mV]')
plt.legend(loc='upper right')
#plt.tight_layout()
#plt.savefig(f'figures/original_izhikevich_model_v_exc_a{a_str}_b{b_str}_c{c_str}_d{d_str}.png', dpi=120)

plt.subplot(2, 1, 2)
plt.plot(v_array[neuron_i_id, :], label=f'inhibitory neuron (id {neuron_i_id})', color='orange')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlim([0, T])
plt.ylabel('$v(t)$ [mV]')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(f'figures/izhikevich_SNN_v_Ne_{name_e}_Ni_{name_i}.png', dpi=120)
plt.show()

plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(u_array[neuron_e_id, :], label=f'recovery variable u(t) (id {neuron_e_id})', color='blue')
plt.title(f'Recovery variable u(t) for typical example neurons')
#plt.title('Recovery Variable u(t) for one example excitatory neuron')
plt.xlabel('Time (ms)')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlim([0, T])
plt.ylabel('$u(t)$ [a.u.]')
plt.legend(loc='upper right')

plt.subplot(2, 1, 2)
plt.plot(u_array[neuron_i_id, :], label=f'recovery variable u(t) (id {neuron_i_id})', color='orange')
#plt.title('Recovery Variable u(t) for one example excitatory neuron')
plt.xlabel('Time (ms)')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.xlim([0, T])
plt.ylabel('$u(t)$ [a.u.]')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(f'figures/izhikevich_SNN_u_Ne_{name_e}_Ni_{name_i}.png', dpi=120)
plt.show()

# %% END