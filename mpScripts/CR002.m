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

% INPUTS   ==============================================================
% source emf  
    emf = 10;
% Resistance values R1   R2    R3   R4   
    R = [0 100 200 400];
% Load resistance R4 = Rload
    RLmin = 10;        
    RLmax = 1000;
    N = 5000;

% CALCULATIONS   ========================================================
   Vm = [emf;0];
   Rm(1,1) = R(1) + R(2 )+ R(3);
   Rm(1,2) = -R(3);
   Rm(2,1) =  R(3);
  
   PLoad = zeros(N,1);
   RLoad = linspace(RLmin,RLmax,N);
   
for c = 1 : N
  R(4) = RLoad(c);
  Rm(2,2) = -(R(3) + R(4));
  Im = Rm \ Vm;
  IR(1) = Im(1);
  IR(2) = Im(1);
  IR(3) = Im(1)- Im(2);
  IR(4) = Im(2);
  V = IR .* R;
  P = IR.^2 .* R;
  PLoad(c) = P(4);
end


% DISPLAY RESULTS   ======================================================
   PLmax = max(PLoad)
   RLmax = RLoad(PLoad == PLmax)
   RLtheory = R(2)*R(3)/(R(2)+R(3))
   
   
% GRAPHICS   =============================================================
figure(1)
   pos = [0.1 0.1 0.25 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = RLoad; yP = 1e3*PLoad;
   plot(xP,yP,'lineWidth',2);

   xlabel('R_{load}  [ \Omega ]');
   ylabel('P_{load}  [ mW ]');
   grid on
   set(gca,'Fontsize',14);
   
