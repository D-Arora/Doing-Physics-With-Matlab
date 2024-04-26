# -*- coding: utf-8 -*-
"""

qmSpectrum.py        April 2024

QUANTUM MECHANICS
Electromagetic visible spectrum


Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH MATLAB 
    https://d-arora.github.io/Doing-Physics-With-Matlab/
    
Reference page for documentation and notes
    https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm006.pdf


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

tStart = time.time()


def colour(wL):
    thiscolor = [0,0,0]
   # wL    = wL*1e+9    # Convert to nm.

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

#%%
# Make visbile spectrum
N = 1501
x = linspace(380,770,N)
y = ones(N)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (6,2)

fig1, axes = plt.subplots(nrows=1, ncols=1)
axes.set_xlabel('$\lambda$  [ nm ]',fontsize = 12)

for c in range(N-2):
    wL = x[c] 
    col = colour(wL)
    xP = x[c:c+1]; yP = y[c:c+1]
    axes.fill_between(xP, yP,color = col)

axes.set_yticks([])
fig1.tight_layout()

fig1.savefig('a1.png')

#%%
tExe = time.time() - tStart
print('  ')
print(r'Execution time:  %2.2f' %tExe)



