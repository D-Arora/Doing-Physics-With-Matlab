% CN06.m

% AC CIRCUIT ANAYSIS MADE SIMPLE
% Solving standard textbook problems with ease
% Each problem is solved within a CELl

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
 
% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180219



%%  CELL 1
clear all
close all
clc

  R = 30
  L = 14e-3
  C = 1e-6
  VS = 20
  f = 1590

  w = 2*pi*f
  ZR = R
  ZC = -1j/(w*C)
  ZL = 1j*w*L
  Z = ZR + ZC + ZL
  
  iS = VS/Z
  IS = abs(iS)
  theta = rad2deg(angle(iS))
  
%% CELL 2
close all
clear all
clc

  R = 300
  L = 2
  VS = 250
  f = 100/pi
  
  w = 2*pi*f
  ZR = R
  ZL = 1j*w*L
  Z = ZR + ZL
  iS = VS/Z
  vS = VS
  vR = iS * ZR
  vL = iS * ZL
  
  IS = abs(iS)
   theta = rad2deg(angle(iS))
   
  VR = abs(vR)
  phiR = rad2deg(angle(vR))
  
  VL = abs(vL)
  phiL = rad2deg(angle(vL))
  
  P = real(iS) * VS
  
  
  figure(1)
    xP = [0 real(vS)]; yP = [0 imag(vS)];
    plot(xP,yP,'b','linewidth',2);
    hold on
    xP = [0 real(vR)]; yP = [0 imag(vR)];
    plot(xP,yP,'r','linewidth',2);
    xP = [0 real(vL)]; yP = [0 imag(vL)];
    plot(xP,yP,'m','linewidth',2);
    grid on
    xlabel('Re'); ylabel('Im');
    legend('v_S','v_R','v_L');
    text(160,-20,'v_S = v_R + v_L','fontsize',12);
    set(gca,'fontsize',12)
  
  
  %% CELL 3
  close all
clear all
clc

  R = 1000
  L1 = 10e-3
  L2 = 250e-3
  C=0.1e-6
  VS = 10
  f = 1000
  
  w = 2*pi*f
  ZR = R
  ZL1 = 1j*w*L1
  ZL2 = 1j*w*L2
  ZC = -1j/(w*C)
  ZP = 1/(1/ZL2 + 1/ZR)
  Z = ZL1 + ZC + ZP
  iS = VS/Z
  IS = abs(iS)
  phi_S = rad2deg(angle(iS))
  vP = iS * ZP
  
  iR = vP / ZR
  IR = abs(iR)
  theta_R = rad2deg(angle(iR))
    
  iL2 = vP / ZL2
  IL2 = abs(iL2)
  theta_L2 = rad2deg(angle(iL2))
  
  figure(1)
    xP = [0 real(iS)]*1e3; yP = [0 imag(iS)]*1e3;
    plot(xP,yP,'b','linewidth',2);
    hold on
    xP = [0 real(iR)]*1e3; yP = [0 imag(iR)]*1e3;
    plot(xP,yP,'r','linewidth',2);
    xP = [0 real(iL2)]*1e3; yP = [0 imag(iL2)]*1e3;
    plot(xP,yP,'m','linewidth',2);
    grid on
    xlabel('Re  [mA]'); ylabel('Im  [mA]');
    legend('i_S','i_R','v_{L2}','location','east');
    text(160,-20,'v_S = v_R + v_L','fontsize',12);
    set(gca,'fontsize',12)
  
 %% CELL 4
 
 clear all
 close all
 clc
 
   vS = 200
   Z1 = 2+1j*4
   Z2 = 2+1j*1
   Z3 = 1-1j*3
   Z4 = 3 + 1j*2
   Z5 = 1/(1/Z1 + 1/Z2 + 1/Z3)
   Z = Z4 + Z5
 
   iS = vS/Z
   IS = abs(iS)
   theta = rad2deg(angle(iS))
   
   p = iS * vS
   thetaP = rad2deg(angle(p))
   
   powerFactor = cosd(thetaP)
   
 