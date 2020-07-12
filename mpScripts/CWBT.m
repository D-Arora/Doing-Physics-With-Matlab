% CR002.m

% Modelling a simple voltage divider circuit
% Maximum energy transfer from source to load

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm
% 171220




clear all
close all
clc

% INPUTS   ===============================================================
   emf = 10;
   R = [100 100 100 100 1e5];

% CALCULATIONS   ======================================================== 
  Vm = [emf;0;0];
  
% Transducer resistance R4 = RTransducer
    RTmin = 99;        
    RTmax = 101;
    N = 5000;
    RT = linspace(RTmin,RTmax,N);
    VT = zeros(N,1);
    
% CALCULATIONS   ========================================================
for c = 1 : N
  R(4) = RT(c);
  Vm = [emf;0;0];
  Rm(1,1) = R(1) + R(2);
  Rm(1,2) = -R(1);
  Rm(1,3) = -R(2);
  
  Rm(2,1) =  R(1);
  Rm(2,2) = -(R(1) + R(3) + R(5));
  Rm(2,3) =  R(5);
  
  Rm(3,1) =  R(2);
  Rm(3,2) =  R(5);
  Rm(3,3) =  -(R(5) + R(4) + R(2));
  
  
  Im = Rm \ Vm;
  
  
  IR(1) = Im(1) - Im(2);
  IR(2) = Im(1) - Im(3);
  IR(3) = Im(2);
  IR(4) = Im(3);
  IR(5) = Im(3) - Im(2);
  
  V = IR .* R;
  VT(c) = V(5); 
end


% DISPLAY RESULTS   ======================================================
  
   
   
% GRAPHICS   =============================================================
figure(1)
   pos = [0.1 0.1 0.25 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = RT; yP = 1e3.*VT;
   plot(xP,yP,'lineWidth',2);

   ylabel('Output Voltage V_5  [ mV ]');
   xlabel('Transducer   R_4  [ \Omega ]');
   grid on
   set(gca,'Fontsize',14);
   
