# -*- coding: utf-8 -*-

"""
Nov 2026         simFM004A.py

FOURIER METHODS - FOURIER TRANSFORM  
  ECG RECORDING

Ian Cooper
 https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/simFM004.pdf

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time
from scipy.integrate import solve_ivp, odeint, simpson
from numpy import pi, linspace, zeros, array, ones, sqrt, cos, sin, exp, real, imag
import scipy.io
from scipy.signal import find_peaks
from scipy.io import wavfile
import os

tStart = time.time()
plt.close('all')

#%%
# Load audio file
def loadFile(name):
    sampleRate, data = wavfile.read(name)
    dt = 1/sampleRate
    return [sampleRate, data, dt]

# Enter 1,2,3,4,5, or 6 to select the sound recording
flag = 6

if flag == 1:   # Pure tone at 440 Hz
   name = 'audio440.wav'               # name of audio file
   sampleRate, data, dt = loadFile(name)   
   t2 = 0.01; t1 = 0    # time span 
   N = int(t2/dt)
   t = linspace(t1,t2,N)
   q = int(100+N); R = np.arange(100,q,1)
   h = data[R]
   h = h/max(h)

   peaks, _ = find_peaks(h)
   peak_times = t[peaks]
   periods = np.diff(peak_times)
   Tavg = np.mean(periods)          # [s]
   favg = 1/Tavg                    # [Hz]
   q = Tavg*1e3; print('period T = %0.2f  ms' % q)
   print('frequency f = %0.0f  Hz' % favg)

if flag == 2:    # beats
    name = 'wav_S1000_1008.wav'               # name of audio file
    sampleRate, data, dt = loadFile(name)   
    t2 = 1; t1 = 0    # time span 
    N = int(t2/dt)
    t = linspace(t1,t2,N)
    q = int(100+N); R = np.arange(100,q,1)
    h = data[R]
    h = h/max(h)
    
    peaks, _ = find_peaks(h)
    peak_times = t[peaks]
    periods = np.diff(peak_times)
    Tavg = np.mean(periods)          # [s]
    favg = 1/Tavg                    # [Hz]
    q = Tavg*1e3; print('period T = %0.2f  ms' % q)
    print('frequency f = %0.0f  Hz' % favg)
 
if flag == 3:    # Guitar  220 Hz
    name = 'audioGuitar1.wav'               # name of audio file
    sampleRate, data, dt = loadFile(name)   
    t2 = 3; t1 = 0    # time span 
    N = int(t2/dt)
    t = linspace(t1,t2,N)
    q = int(100+N); R = np.arange(100,q,1)
    h = data[R,0]
    h = h/max(h)
    f1 = -800; f2 = 800
    
if flag == 4:    # Clarinet  220 Hz
    name = 'audioClarinet1.wav'               # name of audio file
    sampleRate, data, dt = loadFile(name)   
    t2 = 3; t1 = 0    # time span 
    N = int(t2/dt)
    t = linspace(t1,t2,N)
    q = int(100+N); R = np.arange(100,q,1)
    h = data[R,0]
    h = h/max(h)
    f1 = -1200; f2 = 1200
    
if flag == 5:    # Voice
    name = 'audioVoice1.wav'               # name of audio file
    sampleRate, data, dt = loadFile(name)   
    t2 = 3; t1 = 0    # time span 
    N = int(t2/dt)
    t = linspace(t1,t2,N)
    q = int(100+N); R = np.arange(100,q,1)
    h = data[R,0]
    h = h/max(h)
    f1 = -1200; f2 = 1200
    
if flag == 6:    # Train Whistle
    name = 'audioTrain.wav'               # name of audio file
    sampleRate, data, dt = loadFile(name)   
    t2 = 1.5; t1 = 0    # time span 
    N = int(t2/dt)
    t = linspace(t1,t2,N)
    q = int(100+N); R = np.arange(100,q,1)
    h = data[R]
    h = h/max(h)
    f1 = -1800; f2 = 1800    

#%%  FREQUENCY DOMAIN
Nf = 3999
f = linspace(f1,f2,Nf)
p  = 1j*2*pi
HR = zeros(Nf); HI = HR; H = HR+HI*1j
# Fourier transform
for c in range(Nf):
    
     g = h*exp(p*f[c]*t)
     H[c] = simpson(real(g),t) + 1j*simpson(imag(g),t)
# one-sided power spectral density 
psd = 2*abs(H)**2 

# Numerical value of peaks
if flag == 3:
   peaks, _ = find_peaks(abs(H)/max(abs(H)),height = 0.25)              
   print('frequency peaks',f[peaks])


#%% PLOT 1: h(t)     signal
txtX = 't [ ms ]'; txtY = 'h(t) '
plt.rcParams['font.size'] = 12
fig1, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.8))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
xP = t*1e3;
if flag == 3: xP = t
yP = h; ax.plot(xP,yP,'b',lw = 2)
fig1.tight_layout()
fig1.savefig('a1.png')

#%% PLOT 2: H(t)   Fourier transform
txtX = 'f [ Hz ]'; txtY = '| H(f) '
plt.rcParams['font.size'] = 12
fig2, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.8))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
xP = f; yP = abs(H)/max(abs(H)); ax.plot(xP,yP,'b',lw = 2)
fig2.tight_layout()
fig2.savefig('a2.png')

#%% PLOT 3: Fourier transform   psd
txtX = 'f [ Hz ]'; txtY = 'psd '
plt.rcParams['font.size'] = 12
fig3, ax = plt.subplots(nrows=1, ncols=1,figsize = (5,2.8))
ax.set_xlabel(txtX)
ax.set_ylabel(txtY)
ax.grid()
xP = f[f>0]; yP = psd[f>0]/max(psd[f>0]); ax.plot(xP,yP,'b',lw = 2)
fig3.tight_layout()
fig3.savefig('a3.png')

#%%
tExe = time.time() - tStart
print('  ')
print('Execution time = %0.2f s' %tExe)


#%%   Play sound file
#name = 'audioVoice1.wav'
file = name
os.system(file)