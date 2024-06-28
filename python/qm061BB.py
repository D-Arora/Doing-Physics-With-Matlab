# -*- coding: utf-8 -*-
"""
qm061B.py    June 2024

Ian Cooper 
        matlabvisualphysics@gmail.com

Website
       https://d-arora.github.io/Doing-Physics-With-Matlab/
Documentation
       https://d-arora.github.io/Doing-Physics-With-Matlab/pyDocs/qm061.pdf

QUANTUM MECHANICS     
   SPHERICAL HARMONICS   /   ASSOCIATED LEGENDRE FUNCTIONS
   
Origin source of Code
https://scipython.com/blog/visualizing-the-real-forms-of-the-spherical-harmonics/

https://people.csail.mit.edu/sparis/sh/index.php?img=64
   
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# The following import configures Matplotlib for 3D plotting.
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm


#l, m = 5, 2
l = 3
m = 1
plt.rc('text', usetex=True)

# Grids of polar and azimuthal angles
theta = np.linspace(0, np.pi, 220)
phi = np.linspace(0, 2*np.pi, 220)
# Create a 2-D meshgrid of (theta, phi) angles.
theta, phi = np.meshgrid(theta, phi)
# Calculate the Cartesian coordinates of each point in the mesh.
xyz = np.array([np.sin(theta) * np.sin(phi),
                np.sin(theta) * np.cos(phi),
                np.cos(theta)])

def plot_Y(ax, el, m):
    """Plot the spherical harmonic of degree el and order m on Axes ax."""

    Y = sph_harm(abs(m), el, phi, theta)

    # Linear combination of Y_l,m and Y_l,-m to create the real form.
    if m < 0:
        Y = np.sqrt(2) * (-1)**m * Y.imag
    elif m > 0:
        Y = np.sqrt(2) * (-1)**m * Y.real
    Yx, Yy, Yz = np.abs(Y) * xyz

    # Colour the plotted surface according to the sign of Y.
    #cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('PRGn'))
    cmap = plt.cm.ScalarMappable(cmap=plt.get_cmap('jet'))
    cmap.set_clim(-0.5, 0.5)

    ax.plot_surface(Yx, Yy, Yz,
                    facecolors=cmap.to_rgba(Y.real),
                    rstride=2, cstride=2)

    # Draw a set of x, y, z axes for reference.
    ax_lim = 0.8
    ax.plot([-ax_lim, ax_lim], [0,0], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [-ax_lim, ax_lim], [0,0], c='0.5', lw=1, zorder=10)
    ax.plot([0,0], [0,0], [-ax_lim, ax_lim], c='0.5', lw=1, zorder=10)
    # Set the Axes limits and title, turn off the Axes frame.
    #ax.set_title(r'$Y_{{{},{}}}$'.format(el, m))
    
    ax_lim = 0.4
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis('off')
    ax.set_xlabel(r'mmmmm')

#plt.rcParams["figure.figsize"] = (2,1.5)
#fig = plt.figure(figsize=plt.figaspect(1.0))
figsize_px, DPI = 400, 100
figsize_in = figsize_px / DPI
fig = plt.figure(figsize=(figsize_in, figsize_in), dpi=DPI)



ax = fig.add_subplot(projection='3d')

plot_Y(ax, l, m)


fig.tight_layout() 

fig.savefig('a1.png',bbox_inches='tight',dpi=800,
            pad_inches=0.1, transparent=True)
  
#plt.show()

#To plot a family of these functions, try:

# el_max = 3
# figsize_px, DPI = 800, 100
# figsize_in = figsize_px / DPI
# fig = plt.figure(figsize=(figsize_in, figsize_in), dpi=DPI)
# spec = gridspec.GridSpec(ncols=2*el_max+1, nrows=el_max+1, figure=fig)
# for el in range(el_max+1):
#     for m_el in range(-el, el+1):
#         print(el, m_el)
#         ax = fig.add_subplot(spec[el, m_el+el_max], projection='3d')
#         plot_Y(ax, el, m_el)
# plt.tight_layout()
#plt.savefig('a_sph_harm.png')
# plt.show()