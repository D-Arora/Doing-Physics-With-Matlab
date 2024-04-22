# -*- coding: utf-8 -*-
"""

qmSpectrum.py        April 2024

QUANTUM MECHANICS
Blackbody readiation: STARS


Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH MATLAB 
    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qmSun.pdf


Function colour:
    Return the color appropriate to the supplied wavelength.
    Is it assumed the supplied wL is within the range 380-780 nm.
    Smaller or higher values are set notionally to the extreme values. 
    wL input is in nanometre [nm].

"""

#%%  LIBRARIES
from matplotlib import pyplot as plt 
import numpy as np 
from matplotlib.animation import FuncAnimation, PillowWriter  
import time
from numpy import pi, sin, cos, linspace, zeros, ones, exp, sqrt, diag
from scipy.integrate import odeint, quad, dblquad, simps
from scipy.signal import find_peaks

tStart = time.time()





def colour(wL):
    thiscolor = [0,0,0]
    wL    = wL*1e+9    # Convert to nm.

    if wL<380: 
        thiscolor = [1,0,1]

    if (wL>=380)&(wL<440):
        thiscolor = [(440-wL)/(440-380),0,1]

    if (wL>=440)&(wL<490):
        thiscolor = [0,(wL-440)/(490-440),1]

    if (wL>=490)&(wL<510):
        thiscolor = [0,1,(510-wL)/(510-490)]

    if (wL>=510)&(wL<580):
        thiscolor = [(wL-510)/(580-510),1,0]

    if (wL>=580)&(wL<645):
        thiscolor = [1,(645-wL)/(645-580),0]

    if (wL>=645):
        thiscolor = [1,0,0]

#  The intensities fall off near limits of vision

    if wL>700:
       thiscolor = np.array(thiscolor) * (0.3 + 0.7*(780-wL)/(780-700))

    if wL<420:
       thiscolor = np.array(thiscolor) * (0.3 + 0.7*(wL-380)/(420-380))
    
    if thiscolor[0] < 0:
       thiscolor[0] = 0
    if thiscolor[1] < 0:
       thiscolor[1] = 0
    if thiscolor[2] < 0:
       thiscolor[2] = 0   
    return thiscolor

def black(w):
    u1 = K2/w**5; u2 = np.exp(K1/w); u3 = 1/(u2 -1); R = u1*u3
    return R


#%%  Inputs 
# Wavelength at peak of Sun's spectral exitance in wavelenght  5.0225e-7 m
wL_peak = 5.0225e-7          
# Number of data points for wavelength range
num = 801                     

# Constants
c = 2.99792458e8              # speed of light
h = 6.62608e-34               # Planck constant
kB = 1.38066e-23              # Boltzmann constant
sigma = 5.6696e-8             # Stefan constant  
b_wL = 2.898e-3               # Wien constant: wavelength
b_f = 2.82*kB/h               # Wien constant: frequency
R_sun = 6.93e8                # Radius: Sun
R_SE = 1.496e11               # Radius: Sun-Earth
R_E = 6.374e6                 # Radius: Earth
I_0 = 1360                    # Solar constant 
albedo = 0.3                  # Albedo: reflectivity of the Earth's surface
e = 1.602e-19                 # Electron charge

#%%  CALCULATE THE UNKNOWN VARIABLES: Theoretical calculations(SI uniuts)

# Sun's surface temperature  [K]
T = b_wL/wL_peak  
#T = 4000            
# Spectral exitance in wavelength  [W.m-2.m-1]
R_wL = zeros(num)    
# Spectral exitance in frequency  [W.m-2.s-1]  
R_f = zeros(num)            

# Frequency for peak in spectral exitance [Hz]
f_peak = b_f * T    
# Surface area of Sun and Earth [m2]
A_sun = 4*pi*R_sun**2          
A_E = 4*pi*R_E**2              
# Total power output of Sun
P_sun = A_sun * sigma * T**4   
 
#%%  BLACKBODY SPECTRUM  
# R SPECTRAL EXITANCE in wavelength  [W.m-2.m-1] 
K1 = (h*c)/(kB*T)              # constants to simply calculation
K2 = (2*np.pi*h*c**2)

#%% ALL WAVELENGTHS
wL1 = wL_peak/10               # min for wavelength range  lambda1
wL2 = 10*wL_peak               # max for wavelength range  lambda2
wL = linspace(wL1,wL2,num)     # wavelength 
R_wL = black(wL) 
# Wavelength at Peak using logical operations    
wL_peak_graph = wL[np.where(R_wL == max(R_wL))]

# SPECTRAL EXITANCE in visiable wavelength range [W.m-2.m-1] 
wL1_vis = 400e-9;  wL2_vis = 700e-9; num_wL = 95
wL_vis = linspace(wL1_vis,wL2_vis,num_wL) 
R_wL_vis = black(wL_vis) 

# SPECTRAL EXITANCE in infrared wavelength range [W.m-2.m-1] 
wL1_IR = 700e-9;  wL2_IR = wL2
wL_IR = linspace(wL1_IR,wL2_IR,num) 
R_wL_IR = black(wL_IR)  

# SPECTRAL EXITANCE in ultraviolet wavelength range [W.m-2.m-1] 
wL1_UV = wL1;  wL2_UV = 400e-9
wL_UV = linspace(wL1_UV,wL2_UV,num)    
R_wL_UV = black(wL_UV)  

# Area under curves: power output of Sun  [W]
P_total = A_sun * simps(R_wL,wL) 
P_vis   = A_sun * simps(R_wL_vis,wL_vis)
P_IR    = A_sun * simps(R_wL_IR,wL_IR)
P_UV    = A_sun * simps(R_wL_UV,wL_UV)

# Percentage radiation in visible, IR and UV
E_vis = 100 * P_vis / P_total
E_IR  = 100 * P_IR / P_total
E_UV = 100 * P_UV / P_total
  
#%%
# SPECTRAL EXITANCE in frequency  [W.m-2.s-1]  
f1 = f_peak/20             # min for frequency range
f2 = 5*f_peak              # max for frequency range
f = linspace(f1,f2,num)    # frequency
K3 = h/(kB*T)              # constants to simply calcuation
K4 = (2*pi*h/c**2);
R_f = (K4 * f**3) / (np.exp(K3*f)-1) 
  
# Frequency at Peak using logical operations
f_peak_graph = f[R_f == max(R_f)]   

# Area under curve   [W]
P_f = A_sun * simps(R_f,f) 


# SUN-EARTH 
# Solar constant: Intensity of radiation at top of Earth's atmosphere
I_E = P_total/(4*np.pi*R_SE**2)  # [W.m-2]   
# Power absorbed by Earth from Sun  [W]
P_E = (1-albedo) * I_E * np.pi* R_E**2   
# Temperature of the Earth  [K]
T_E = ((1-albedo)*I_E/(4*sigma))**0.25 
Tc = T_E - 273

#%%  CONSOLE OUTPUT
print('')
print('Sun: temperature of photosphere, T_S = %2.0f K ' %T)
print('Peak in Solar Spectrum')
print('  Theory: Wavelength at peak in spectral exitance, wL_peak = %2.2e m' %wL_peak)
print('  Graph: Wavelength at peak in spectral exitance,  wL_peak = %2.2e m' %wL_peak_graph)
print('  Theory: Frequency at peak in spectral exitance, f_peak = %2.2e Hz' %f_peak)
print('  Graph:  Frequency at peak in spectral exitance, f_peak = %2.2e Hz' %f_peak_graph)
print('Total Solar Power Output')
print('  P_Stefan_Boltzmann = %2.2e  W' % P_sun)
print('  P_wL               = %2.2e  W' % P_total)
print('  P_f                = %2.2e  W' % P_f)
print('IR / visible / UV')
print('  P_IR  = %2.2e W   ' %P_IR +'percentage %2.2f' %E_IR)
print('  P_vis = %2.2e W   ' %P_vis +'percentage %2.2f' %E_vis)
print('  UV    = %2.2e W   ' %P_UV +'percentage %2.2f' %E_UV)
print('  ')
print('Sun - Earth')   
print('Theory: Solar constant I_O   = 1.360e+03  W/m^2') 
print('Computed: Solar constant I_E = %2.2e w/m^2' %I_E) 
print('Surface temperature of the Earth, T_E  = %2.0f K' %T_E + ' = %2.0f deg C ' %Tc)
  
#%%
# GRAPHICS

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,4)

fig1, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('$\lambda$  [ m ]',fontsize = 12)
#axes.set_ylabel('$R_{\lambda}$  [ W.m$^{-2}.m$^{-1}$ ]',fontsize = 12)
axes.set_ylabel('R$_{\lambda}$  [ W.m$^{-2}.m^{-1}$ ] ',fontsize = 12)
#axes.set_ylim([0,10e-7])
#axes.set_xlim([0,3e-6])
for s in range(num-2):
    w = wL[s] 
    col = colour(w)
    xP = wL[s:s+1]; yP = R_wL[s:s+1]
    axes.fill_between(xP, yP,color = col)
    fig1.tight_layout()

fig1.savefig('a1.png')


#%%
# GRAPHICS

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,4)

fig1, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('$f $  [ Hz ]',fontsize = 12)
#axes.set_ylabel('$R_{\lambda}$  [ W.m$^{-2}.m$^{-1}$ ]',fontsize = 12)
axes.set_ylabel('R$_f$  [ W.m$^{-2}.m^{-1}$ ] ',fontsize = 12)
#axes.set_ylim([0,10e-7])
#axes.set_xlim([0,5e-6])
for s in range(num-2):
    w = c/f[s] 
    col = colour(w)
    xP = f[s:s+1]; yP = R_f[s:s+1]
    axes.fill_between(xP, yP,color = col)
    fig1.tight_layout()

fig1.savefig('a2.png')

#%%
tExe = time.time() - tStart
print('  ')
print(r'Execution time:  %2.2f' %tExe)



