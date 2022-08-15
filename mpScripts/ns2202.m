% ns2202.m



clear
close all
clc


%  INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% External current density   [A/cm^2]
%  Jext = 1e-4;
% Specific conductances     [S/cm^2]
  gNa = 0e-3;
  gK  = 36e-3;
  gL = 3e-3;
% Resting potneital  [V]
  Vrest = -0e-3;
% Time interval  [ms]
  tMax = 15e-3; 
  num = 10001;

% SETUP  =========================================================
% Leak conductance  [S/cm^2]
gL =  0.3e-3;
% Specific capacitance  [F/cm^2]
  cm = 1e-4;
% Time
  t = linspace(0,tMax,num);
  dt = t(2) - t(1);
% Reversal potentials  [V]
  ENa =  50e-3;
  EK  = -77e-3;
  EL  = -55e-3;

V = zeros(num,1); JL = V; JK = V; JNa = V; Jm = V;
V(1) = Vrest;

 Jext = 1e-2.*ones(num,1);

% Solve ODE
  for n = 1 : num-1
   JL(n)  = gL* (V(n) - EL);
   JK(n)  = gK* (V(n) - EK);
   JNa(n) = gNa*(V(n) - ENa);
   Jm(n)  = JL(n) + JK(n) + JNa(n); 
   V(n+1) = V(n) + (dt/cm)*(Jext(n) - Jm(n));
  end


  figure(1)
      plot(t.*1e3,V.*1e3)



  %%
clear
close all
clc

V1 = -10
V2 = 20
R1 = 10
R2 = 10

I = (V2 + V1)/(R1+R2)

Vm2 = V2 - I*R2
Vm1 = V1 - I*R1
