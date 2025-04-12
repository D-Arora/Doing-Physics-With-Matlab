# -*- coding: utf-8 -*-
"""
emGB01.py          mar 2025

COMPUTATIONAL OPTICS
   GAUSSIAN BEAMS  (PARAXIAL REGIME)


Ian Cooper
    https://d-arora.github.io/Doing-Physics-With-Matlab/

Documentation
    https://d-arora.github.io/Doing-Physics-With-Matlab/emGB01.pdf
"""

#%%
import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax, array, sqrt, arctan 
import matplotlib.pyplot as plt
import matplotlib.pyplot as mpl
import time
from mpl_toolkits.mplot3d import Axes3D


#%% INPUTS
# Max power transmitted in beam  [W]
P0 = 1e-3
# beam waist  [m]
w0 = 0.5e-3
# wavelength  [m]    [HeNe laser  wL = 632.8e-9 m]
L1 = 632.8e-9; col1 = [1,0,0]
L2 = 480e-9;   col2 = [0,0,1]

# Length of Z domain in multiples of zR :  z = nZ * zR
nZ = 5
# Number of grid points in calculations  
N = 501
# Z position of XY plane for radial irradiance plots:  zPR = zP / zR
#   [e,g,  zPR = 0, 1, 2, ... ]
zPR = 0;

# Speed of light
c = 299792458
# Permittivity of free space
eps0 = 8.85418782e-12

eps = 1e-16

#%%   CALCULATIONS
# Wave number  [1/m]
k1 = 2*pi/L1; k2 = 2*pi/L2
# Angular velocity
omega1 = c*k1; omega2 = c*k2
# Rayleigh range  [m]
zR1 = pi*w0**2/L1; zR2 = pi*w0**2/L2

# # Electric field amplitude  [V/m]
E0 = sqrt(4*P0/(pi*c*eps0*w0**2))

      
# # Max beam irradiance  [W/m^2]
# Smax = 2*P0/(pi*w0**2)
A = (c*eps0*E0**2/2)      # same as Smax
# z domain
z = linspace(0,5,N)

# Beam divergence angle [rad]  [deg]
thetaR1 = w0/zR1; thetaR2 = w0/zR2;
theta1 = 180*thetaR1/pi; theta2 = 180*thetaR2/pi;


#%% GRAPHICS
# 1   Beam spot w  [m]
w1 = w0 * sqrt(1+(z/zR1)**2); w2 = w0 * sqrt(1+(z/zR2)**2)

plt.rcParams["figure.figsize"] = (6,3)
fig1, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('z  [m]',fontsize = 12)
ax.set_ylabel('beam spot  w  [ mm ] ',fontsize = 12)
ax.text(0.2,2.3,'633 nm', color = col1)
ax.text(0.2,2.1,'480 nm', color = col2)
ax.text(4,2.2,'$\Theta$ = %0.3f deg ' %theta1, color = col1)
ax.text(4,1.1,'$\Theta$ = %0.3f deg ' %theta2, color = col2)
ax.text(1.8,2.2,'$z_R$ = %0.3f m ' %zR2, color = col2)
ax.text(0.1,1.2,'$z_R$ = %0.3f m ' %zR1, color = col1)
ax.grid()
xP = z; yP = w1*1e3; ax.plot(xP,yP,color = col1,lw = 2) 
yP = w2*1e3; ax.plot(xP,yP,color = col2,lw = 2) 
xP = [zR1,zR1]; yP = [0,2.5]
ax.plot(xP,yP,'r',lw = 1)
xP = [zR2,zR2]; yP = [0,2.5]
ax.plot(xP,yP,'b',lw = 1)
 
xP = [0,5]; yP = [sqrt(2)*w0*1e3,sqrt(2)*w0*1e3]; ax.plot(xP,yP,'k',lw = 2) 

fig1.tight_layout()

# 2   Radius of curvature
R1 = (z + eps) + zR1**2/(z + eps); R2 = (z + eps) + zR2**2/(z + eps);
plt.rcParams["figure.figsize"] = (6,3)
fig2, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('z / z$_R$',fontsize = 12)
ax.set_ylabel('R [ m ] ',fontsize = 12)
ax.text(3,1.8,'633 nm', color = col1)

ax.text(3,6,'480 nm', color = col2)
ax.text(2.2,15.5,'$z_R$ = %0.3f m ' %zR2, color = col2)
ax.text(2.2,18,'$z_R$ = %0.3f m ' %zR1, color = col1)
ax.grid()
ax.set_ylim([0,20])
xP = z; yP = R1; ax.plot(xP,yP,'r',lw = 2) 
xP = z; yP = R2; ax.plot(xP,yP,'b',lw = 2) 
xP = [zR1,zR1]; yP = [0,20]; ax.plot(xP,yP,'r',lw = 1) 
xP = [zR2,zR2]; yP = [0,20]; ax.plot(xP,yP,'b',lw = 1) 
fig2.tight_layout()

# 3   Gouy phase
plt.rcParams["figure.figsize"] = (6,3)
fig3, ax = plt.subplots(nrows=1, ncols=1)
Z = linspace(-z,z,2*N)
GP1 = arctan(Z/zR1)/pi; GP2 = arctan(Z/zR2)/pi;
ax.set_xlabel('z  [ m ]',fontsize = 12)
ax.set_ylabel('$\phi / \pi$ ',fontsize = 12)
# # q = wL*1e9; ax.set_title(r'$\lambda$ = %0.0f  nm   ' %q + 'z$_R$ = %0.3f  m' %zR + 
# #              '   $\Theta$ = %0.3f deg ' %theta, fontsize = 12)
ax.grid()
ax.set_xlim([-5,5]); ax.set_ylim([-0.5,0.50])
ax.set_xticks(np.arange(-5,5.5,1))
ax.set_yticks(np.arange(-0.5,0.55,0.1))
ax.text(-4.4,-0.08,'633 nm', color = col1)
ax.text(-4.4,-0.18,'480 nm', color = col2)
xP = Z; yP = GP1; ax.plot(xP,yP,'r',lw = 2) 
xP = Z; yP = GP2; ax.plot(xP,yP,'b',lw = 2) 
fig3.tight_layout()

# 4   Axial radiance
Sz1 = A / (1+(z/zR1)**2)
Sz2 = A / (1+(z/zR2)**2)

plt.rcParams["figure.figsize"] = (6,3)
fig4, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('z / z$_R$',fontsize = 12)
ax.set_ylabel('$S_z$  [ W.m$^{-2}$ ] ',fontsize = 12)
ax.text(2.2,2300,'633 nm', color = col1)
ax.text(2.2,2050,'480 nm', color = col2)
ax.grid()
ax.set_xlim([0,5])
ax.set_xticks(np.arange(0,5.5,1))
ax.set_ylim([0,2700])
xP = z; yP = Sz1; ax.plot(xP,yP,'r',lw = 2) 
xP = z; yP = Sz2; ax.plot(xP,yP,'b',lw = 2) 
xP = [zR1,zR1]; yP = [0,2500]; ax.plot(xP,yP,'r',lw = 1) 
xP = [zR2,zR2]; yP = [0,2500]; ax.plot(xP,yP,'b',lw = 1) 
xP = [0,5]; yP = [max(Sz1)/2,max(Sz1)/2]; ax.plot(xP,yP,'k',lw = 1) 
fig4.tight_layout()

# 5   Radial radiance    [W/m^2]
zP = 5
wP1 = w0 * sqrt(1+(zP/zR1)**2)
wP2 = w0 * sqrt(1+(zP/zR2)**2)
r = linspace(0,10*w0,N)       # r/w

A = (c*eps0*E0**2/2)
Sr1 = A*(w0/wP1)**2*exp(-2*(r/wP1)**2) 
Sr2 = A*(w0/wP2)**2*exp(-2*(r/wP2)**2) 

plt.rcParams["figure.figsize"] = (6,2.5)
fig5, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('r / w$_0$',fontsize = 12)
ax.set_ylabel('$S_r$  [ W.m$^{-2}$ ] ',fontsize = 12)
q = amax(array([Sr1,Sr2]))
ax.text(5.2,0.4*q,'633 nm', color = col1)
ax.text(5.2,0.6*q,'480 nm', color = col2)
ax.text(5.05,0.8*q,'z$_P$ = %0.2f  m' %zP, color = 'k',fontsize = 12)
ax.grid()
xP = r/w0; yP = Sr1; ax.plot(xP,yP,'r',lw = 2) 
xP = r/w0; yP = Sr2; ax.plot(xP,yP,'b',lw = 2) 
fig5.tight_layout()


#%%   6   2D / 3D intensity plots in XY plane   [W/m^2]
zP = 5
rP = 2e-3
flag = 633      # 633 or 480

x = linspace(-rP,rP,N)
[xx, yy] = np.meshgrid(x,x)
r2 = xx**2 + yy**2

if flag == 633:
    Col = 'Reds' 
    zR = zR1
    wP = w0 * sqrt(1+(zP/zR)**2)
if flag == 480:
    Col = 'Blues' 
    zR = zR2
    wP = w0 * sqrt(1+(zP/zR)**2)

K = -2*r2/wP**2
Sxy = A*(w0/wP)**2 * exp(K)

plt.rcParams["figure.figsize"] = (3.5,3.5)
fig6, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('x  [ mm ]',fontsize = 12)
ax.set_ylabel('y  [ mm ]',fontsize = 12)
ax.set_title('max(S$_{xy}$) = %0.0f  W.m$^{-2}$' %amax(Sxy))
xP = xx*1e3; yP = yy*1e3
ax.pcolor(xP,yP,Sxy,cmap = Col)
ax.set_aspect('equal', 'box')
fig6.tight_layout()
 

#%%  7 XY Intensity
plt.rcParams["figure.figsize"] = (3,4)
fig7 = plt.figure()
ax = fig7.add_subplot(1, 1, 1, projection='3d')
ax.set_zticks([])
ax.set_xlabel('x  [ mm ]',fontsize = 10)
ax.set_ylabel('y  [ mm ]',fontsize = 10)
ax.set_title('max(S$_{xy}$) = %0.0f  W.m$^{-2}$' %amax(Sxy),fontsize = 10)
xP = xx*1e3; yP = yy*1e3
ax.set(zlabel=None)
ax.plot_surface(xP,yP,Sxy,cmap = Col)
fig7.tight_layout()


#%% 8 XZ Intensity    [W/m^2]
zP = linspace(0,6,N)
rP = 2e-3*linspace(-1,1,N)
ZP, RP = np.meshgrid(zP,rP)

flag = 480     # 633 or 480


if flag == 633:
    Col = 'Reds'; col = [1,0,0] 
    zR = zR1
    wP = w0 * sqrt(1+(ZP/zR)**2)
if flag == 480:
    Col = 'Blues'; col = [0,0,1] 
    zR = zR2
    wP = w0 * sqrt(1+(ZP/zR)**2)

w = w0*sqrt(1+(zP/zR)**2)
K = -2*RP**2/wP**2
Szr = A*(w0/wP)**2 * exp(K)

plt.rcParams["figure.figsize"] = (6,3)
fig8, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('z  [ m ]',fontsize = 12)
ax.set_ylabel('r  [ mm ]',fontsize = 12)
ax.set_title('max(S$_{zr}$) = %0.0f  W.m$^{-2}$' %amax(Szr))
ax.set_xlim([0,5])
xP = zP; yP = RP*1e3
cf = ax.pcolor(xP,yP,Szr,cmap = Col)

fig8.colorbar(cf, ax=ax)

xP = z; yP = w*1e3;  ax.plot(xP,yP,color = col,lw = 2) 
xP = z; yP = -w*1e3; ax.plot(xP,yP,color = col,lw = 2) 
fig8.tight_layout()


#%%  XY power    red light L = 633 nm
zP = 5
wP = w0*sqrt(1 + (zP/zR1)**2)
rP = linspace(0,10*w0,N)     # rw = rP/wP
q = -2*(rP/wP)**2
P = 100*(1 - exp(q))
Pw = P[wP>rP][-1]

plt.rcParams["figure.figsize"] = (6,3.5)
fig9, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('aperature radius  [ mm ]',fontsize = 12)
ax.set_ylabel('% power thru. aperture',fontsize = 12)
ax.set_title('power thru. beam spot = %0.1f  ' %Pw)
xP =rP*1e3; yP = P
ax.plot(xP,yP,'r',lw = 2)
xP = array([wP,wP])*1e3; yP = [0,Pw]; ax.plot(xP,yP,'k',lw = 1)
xP = array([0,wP])*1e3;  yP = [Pw,Pw]; ax.plot(xP,yP,'k',lw = 1)
ax.grid()
fig9.tight_layout()


#%% Gouy phase
q = linspace(-5,5,N) 
phi = arctan(q)
plt.rcParams["figure.figsize"] = (6,3.5)
fig10, ax = plt.subplots(nrows=1, ncols=1)
ax.set_xlabel('z / z$_R$ ',fontsize = 12)
ax.set_ylabel('$\phi / \pi$',fontsize = 12)
ax.set_xticks(np.arange(-5,6,1))
ax.set_yticks(np.arange(-0.5,0.7,0.25))
xP =q ; yP = phi / pi
ax.plot(xP,yP,'r',lw = 2)
xP = array([1,1]); yP = [-0.5,0.5]; ax.plot(xP,yP,'k',lw = 1)
xP = array([-1,-1]); yP = [-0.5,0.5]; ax.plot(xP,yP,'k',lw = 1)
ax.grid()
fig10.tight_layout()


#%%
# fig1.savefig('a1.png')
# fig2.savefig('a2.png')
# fig3.savefig('a3.png')
# fig4.savefig('a4.png')
# fig5.savefig('a5.png')
# fig6.savefig('a6.png')
# fig7.savefig('a7.png')
# fig8.savefig('a8.png')
# fig9.savefig('a9.png')
# fig10.savefig('a10.png')