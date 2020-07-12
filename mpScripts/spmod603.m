% cemVE04.m
% 22 aug 2017
% Ian Cooper
% School of Physics, University of Sydney
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% calculations for an electron and a proton in a uniform electric field


clear all
close all
clc

% input parameters
   me = 9.11e-31;
   mp = 1.67e-27;
   q = 1.602e-19;
   dV = 1000;
   dd = 20e-3;
   s =  5e-3;
 % Calculations
 E = dV / dd
 F = q*E
 RM = mp/me
 ae = F/me
 ap = F/mp
 te = sqrt(2*s/ae)
 tp = sqrt(2*s/ap)
 ve = ae*te
 vp = ap*tp
 W = F*s
 Ke = 0.5*me*ve^2
 Kp = 0.5*mp*vp^2
 pe = me*ve
 pp = mp*vp
 
 disp('charges in B-field radius of orbit')
  v = 8e5
  B = 0.005
  
  Re = me * v / (q * B)
  Rp = mp * v / (q * B)
  
 
