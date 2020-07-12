% Cac3.m

% ac circuit modeling and analysis
% Solving textbook style problesm on RLC circuits

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180110 



%% CELL #1  Textbook Example 1
clear all
close all
clc

% INPUTS [SI Units] ---------------------------------------------------
R = 100;
C = 10e-6;
L = 100e-3;
f = 60;
VS = 141.4;

% CALCULATIONS  [SI Units] -------------------------------------------
w = 2*pi*f
XC = 1/(w*C)
XL = w*L
Z1 = R;
Z2 = -1i*XC
Z3 = 1i*XL
Z = Z1+Z2+Z3
IS = VS/Z
Ipeak = abs(IS)
Irms = abs(IS)/sqrt(2)
theta = angle(IS)
V1 = IS * Z1
V2 = IS * Z2
V3 = IS * Z3
V23 = V2+V3
V123 = V1+V2+V3
V1rms = abs(V1)/sqrt(2)
V23rms = abs(V23)/sqrt(2)
V123rms = abs(V123)/sqrt(2)





%%   CELL #2    Textbook Example 2
   clear all
   close all
   clc
% INPUTS           SI units
   R = 1000;
   C = 10e-6;
   L = 100e-3;
   f = 60;
   Vin = 100;
 
% CALCULATIONS    SI units 
   w = 2*pi*f;
   ZR = R;
   XC = 1/(w*C);
   XL = w*L;
   ZC = -1j * XC;
   ZL =  1j * XL;
   Z = ZR + ZC +ZL;
   Iin = Vin / Z;
   VR = Iin * ZR;
   VC = Iin * ZC;
   VL = Iin * ZL;
   phiC = angle(VC)/pi;
   phiL = angle(VL)/pi;
   f0 = 1/(2*pi*sqrt(L*C));
   Pin = Vin * conj(Iin);
   PR  = VR  * conj(Iin);
   PC  = VC  * conj(Iin);
   PL  = VL  * conj(Iin);
   
 % DISPLAY RESULTS  actual (real) values
   disp('Inputs  [SI Units]  ');
     fprintf('  R =  %3.2f  \n',R);
     fprintf('  C =  %3.2e  \n',C);
     fprintf('  L =  %3.2e  \n',L);
     fprintf('  f =  %3.2f  \n',f);
     fprintf('  Vin =  %3.2f  \n',Vin);
   disp('Calculations  [SI UNITS]  ')
     fprintf('  XC =  %3.2f  \n',XC);
     fprintf('  XL =  %3.2f  \n',XL);
     fprintf('  Iin =  %3.2f  \n',abs(Iin));
     fprintf('  VR =  %3.2f  \n',abs(VR));
     fprintf('  VC =  %3.2f  \n',abs(VC));
     fprintf('  phiC/pi =  %3.2f  \n',phiC)
     fprintf('  VL =  %3.2f  \n',abs(VL));
     fprintf('  phiL/pi =  %3.2f  \n',phiL);
     fprintf('  VR + VC + VL =  %3.2f  \n',abs(VR+VC+VL));
     fprintf('  VL =  %3.2f  \n',abs(VL));
     fprintf('  phiL/pi =  %3.2f  \n',phiL);
     fprintf('  VR + VC + VL =  %3.2f  \n',abs(VR+VC+VL));
  disp('  ');
     fprintf('  resonance frequency f0 =  %3.2f  \n',f0);
   disp('  ');
     fprintf('  Pin =  %3.2f  \n',abs(Pin));
     fprintf('  PR  =  %3.2f  \n',abs(PR));
     fprintf('  PC  =  %3.2f  \n',abs(PC));
     fprintf('  PL  =  %3.2f  \n',abs(PL));
     fprintf('  PR + PC + PL  =  %3.2f  \n',abs(PR+PC+PL));
     
     
     
%%  CELL #3
 clear all
 close all
 clc
 format shorte
 % INPUTS
   R = 1000;
   C = 10e-6;
   L = 100e-3;
   f = 50;
  Vin = 100;
 
% CALCULATIONS   
   w = 2*pi*f;
   ZR = R;
   ZC = -1j / (w*C);
   ZL = 1j * w*L;
   Zp = 1/(1/ZC + 1/ZL);
   Z = ZR + Zp;
   Iin = Vin / Z;
   IR = Iin;
   VR = IR * ZR;
   Vp = Vin - VR;
   IC = Vp / ZC;
   IL = Vp / ZL;
   f0 = 1/sqrt(L*C);
   
   thetaR = angle(IR);
   thetaC = angle(IC);
   thetaL = angle(IL);
   
   % Display results  actual values (real)
    disp('Inputs  ');
      fprintf('R =  %3.2f  ohms  \n',R);
      fprintf('C =  %3.2e  F  \n',C);
      fprintf('L =  %3.2e  H  \n',L);
      fprintf('f =  %3.2f  Hz  \n',f);
      fprintf('Vin =  %3.2f  V  \n',Vin);
    disp('Calculations  ');
      fprintf('Vin =  %3.2f  V  \n',Vin);
      fprintf('VR =  %3.2f   V \n',abs(VR));
      fprintf('Vp =  %3.2f  V  \n',abs(Vp));
      fprintf('VR + Vp =  %3.2f  V  \n',abs(VR + Vp));
      disp('  ');
      fprintf('IR =  %3.2f  mA  \n',1e3*abs(IR));
      fprintf('IC =  %3.2f  mA  \n',1e3*abs(IC));
      fprintf('IL  =  %3.2f  mA  \n',1e3*abs(IL));
      fprintf('IC + IL =  %3.2f  \n',1e3*abs(IC+IL));
      fprintf('thetaR / pi =  %3.2f  \n',thetaR/pi);
      fprintf('thetaC / pi=  %3.2f  \n',thetaC/pi);
      fprintf('thetaL / pi  =  %3.2f  \n',thetaL/pi);
      fprintf('resonance frequency f0 =  %3.2f Hz  \n',f0);
   
   
   
   
 
 
 
