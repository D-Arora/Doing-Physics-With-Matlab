# -*- coding: utf-8 -*-
"""

qm001.py
April 2024

QUANTUM MECHANICS
Finite Difference Time Development Method
[2D] Schrodinger Equation
Propagation of a [2D] Gausssian Pulse
Arbitrary units are used
Free propgation / Potential Hill / Potential Cliff
Single Slit / Double Slit
The variable flagU (1,2,3,4,5) is used to change the potential energy

Animation can be saved as a gif file (flag
Script - modified version(Kevin Berwick)
        'Computational Physics' by Giordano Nakanishi

Ian Cooper
email: matlabvisualphysics@gmail.com

DOING PHYSICS WITH MATLAB 
    
Reference page for documentation and notes
    

"""

# https://www.tutorialspoint.com/how-to-animate-a-pcolormesh-in-matplotlib

import numpy as np
from numpy import pi, sin, cos, exp, linspace, zeros, amax 
from numpy.linalg import eig
from scipy.integrate import odeint, quad, dblquad, simps
from scipy import pi, sqrt
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import axes3d
import sympy as sym


tStart = time.time()



# Setup [2D] grid  ====================================================
# Position & Time increments / Number of time Steps Constant in S.E.
N = 99
nT = 2

xG = linspace(0,1,N)
xP, yP = np.meshgrid(xG,xG)
dx = xG[2] - xG[1]
dt = 1e-5
f =  dt/(2*dx**2) 
  
#x = np.arange(1,N-2)
#y = x
  
# Potential Energy Function  ==========================================


#   flagU = 4;
U = zeros([N,N])
#   G = round(N/2);
#   switch flagU
#      case 1                           % Free propagation
   
#      case 2
#        U(:,G:N) =  1e3;      % Potential Hill 
#      case 3
#        U(:,G:N) = -1e3;      % Potential Well  
#      case 4                  % Single Slit  
#        U(1:G-5,G-3:G+3)   = 15e3;
#        U(G+5:N,G-3:G+3)   = 15e3;
#      case 5                  % Double Slit  
#        U(1:G-10,G-3:G+3)   = 15e3;
#        U(G+10:N,G-3:G+3)   = 15e3;  
#       U(G-5:G+5,G-3:G+3)   = 15e3;  
#   end
  
 
#[2D] GAUSSIAN PULSE (WAVE PACKET) ==================================
# Initial centre of pulse
x0 = 0.20;  y0 = 0.5
# Initial amplitude of pulse 
A = 10
# Pulse width: sigma squared
s = 5e-3
# Wavenumber
k0 = 50 

# Envelope
psiE = A*exp(-(xP-x0)**2/s)*exp(-(yP-y0)**2/s)
# Plane wave propagation in +X direction
psiP = exp(1j*k0*xP)
# Wavefunction
psi1 = psiE*psiP
# Probability Density  
prob1 = np.conj(psi1)*psi1
# Extract Real and Imaginary parts
R1 = np.real(psi1);  I1 = np.imag(psi1)
prob1 = R1**2 + I1**2
 
plt.pcolor(xP,yP,np.real(prob1))
   
# UPDATE WAVEFUNCTION & GRAPHICS  ====================================
# si1 (current value) psi2 (next value)

plt.rcParams['font.size'] = 12
plt.rcParams["figure.figsize"] = (3.9,3)


fig, ax = plt.subplots(nrows=1, ncols=1)

R2 = zeros([N,N])
I2 = zeros([N,N])
probD = zeros([N,N])

for c in range(1,2):
  #breakpoint()
  for x in range(1,N-1):
      for y in range(1,N-1):
          mR1 = R1[x,y]
          mR2 = -f*(I1[x+1,y]-2*I1[x,y]+I1[x-1,y]+I1[x,y+1] -2*I1[x,y]+I1[x,y-1])
          mR3 = dt*U[x,y]*I1[x,y]
         # R2[x,y] = R1[x,y] - f*(I1[x+1,y]-2*I1[x,y]+I1[x-1,y]+I1[x,y+1] /
         #                -2*I1[x,y]+I1[x,y-1]) + 
          R2[x,y] = mR1 + mR2 + mR3
          
          R1[x,y] = R2[x,y]

          I2[x,y] = I1[x,y] + f*(R1[x+1,y]-2*R1[x,y]+R1[x-1,y]+R1[x,y+1] /
                         -2*R1[x,y]+R1[x,y-1]) - dt*U[x,y]*R1[x,y]
  
          probD[x,y] = R2[x,y]**2 + I1[x,y]*I2[x,y]
          I1[x,y] = I2[x,y]
  
plt.pcolor(xP,yP,np.real(probD))
#  ax.pcolormesh(xP,yP,probD)
  # plt.show()
  #time.sleep(0.5)
  
# subplot(2,1,1) 
# % Graph only updated after a specified number of time steps
# % Probability Density values scaled to shoow very small values 
# if rem(c, 8) == 0
#    surf(x,y, abs(prob2.^0.05))    
#    title('Probability Density (scaled)','fontweight','normal')
#    xlabel('x')
#    ylabel('y')
#    zlabel('ps*psi');
#    axis([0 1 0 1 0 2])
#    view(44,55)
#    axis off
#    light
#    lighting phong
#    camlight('left')
#    shading interp
#   % colorbar
   
# subplot(2,1,2)
#    pcolor(x,y, abs((prob2).^0.05));
#    %title('Probability density','fontweight','normal');
#    xlabel('x')
#    ylabel('y')
#    axis off
#    shading interp
#  % colorbar
#    axis square
#    pause(0.00001)
# end

#    if flagG > 0 
#          frame = getframe(1);
#          im = frame2im(frame);
#          [imind,cm] = rgb2ind(im,256);
#        %  On the first loop, create the file. In subsequent loops, append.
#          if nt == 9
#           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
#          else
#          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
#          end
#          nt = nt+1;
#     end   
# end


# figure(2)
#   set(gcf,'units','normalized');
#   set(gcf,'position',[0.40 0.1 0.25 0.45]);
#   set(gcf,'color','w');
 
# subplot(2,1,1)
#     pcolor(x,y,prob1)
#     shading interp
#     set(gca,'xtick',-100:100)
#     set(gca,'ytick',-100:100)
#     axis square
#     title('Initial Probability Density','fontweight','normal')
#     xlabel('x')
#     ylabel('y')
#     set(gca,'fontsize',10)
    
#  subplot(2,1,2)
#     surf(x,y,prob1)
#     view(34,18)
#     shading interp
#     set(gca,'xtick',-100:100)
#     set(gca,'ytick',-100:100)
#     xlabel('x')
#     ylabel('y')
#     set(gca,'fontsize',10)
 
    
#  figure(3)
#   set(gcf,'units','normalized');
#   set(gcf,'position',[0.7 0.1 0.2 0.30]);
#   set(gcf,'color','w');

#   if flagU == 4 || flagU == 5
#     pcolor(x,y,U)
#   else
#     surf(x,y,U)
#     view(34,18)
#   end
#   shading interp
#   hold on
  
#   set(gca,'xtick',-100:100)
#   set(gca,'ytick',-10000:10000)
#   axis square
#   title('Potential Energy Function','fontweight','normal')
#   xlabel('x')
#   ylabel('y')
#   set(gca,'fontsize',10)
   
    
