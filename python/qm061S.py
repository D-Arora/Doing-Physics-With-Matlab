# -*- coding: utf-8 -*-
"""
qm061C.py    June 2024

Ian Cooper 
        matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm061.pdf

QUANTUM MECHANICS     
   SPHERICAL HARMONICS ON A SPHERE
   
Origin source of Code
    https://scipython.com/book/chapter-8-scipy/examples/visualizing-the-spherical-harmonics/

https://irhum.github.io/blog/spherical-harmonics/index.html
"""

import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import sph_harm


#%% INPUTS
L = 2         # orbital angular momentum number
mL = 2         # azimuthal (magentic) quantum number


phi = np.linspace(0, np.pi, 199)
theta = np.linspace(0, 2*np.pi, 199)
phi, theta = np.meshgrid(phi, theta)

# The Cartesian coordinates of the unit sphere
x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)


# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
fcolors = sph_harm(mL, L, theta, phi).real
fmax, fmin = fcolors.max(), fcolors.min()
fcolors = (fcolors - fmin)/(fmax - fmin)

# Set the aspect ratio to 1 so our sphere looks spherical
fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')
ax.set_title('l = %2.0f' %L + '   $m_l = %2.0f$' %mL)
ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
# Turn off the axis planes



ax.set_axis_off()
ax.set_aspect('equal')

fig.savefig('a1.png',bbox_inches='tight',dpi=800,
             pad_inches=0.1, transparent=True)



#%%

