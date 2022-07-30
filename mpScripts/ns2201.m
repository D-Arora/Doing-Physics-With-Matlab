% ns220104.m

% MODELLING A NEURON  RC circuit

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220723 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%  https://d-arora.github.io/Doing-Physics-With-Matlab/

close all
clear
clc

tic

% INPUTS  =============================================================
% default values and units  [ ]

% capacitance [ C = 1e-10  F]
   C = 1e-10; 
% resistance [R  = 1e8]  
   R = 1e8;
% Initial conditions [V = 0 mV]
   V0 = -20e-3;
% Reversal potential [mV]
  E1 = -75e-3;
  E2 = 55e-3;


% Time span  [0 20 s]
  N = 51;
  tm = linspace(0,200,N);    % time interval [ms]
  t = tm./1e3;                 % time interval [s]

  tSpan = t;

% Relative tolerance for ODE solver  [1e-6]
  RelTol = 1e-6;

% External current   [ms  A]
  tON = 20; tOFF = 80; Imax = 1e-9;
  Iext = zeros(N,1);
  z1 = find(tm > tON,1);  z2 = find(tm > tOFF,1);
  Iext(z1:end) = Imax; 
  Iext(z2:end) = 0;   
  
% CALCULATIONS  ======================================================= 

% time constant [s] % [ms]
  tau = R*C;
  taum = tau*1e3;

  K(1) = C; K(2) = R;
  K(3) = tON; K(4) = tOFF; K(5) = Imax;
  K(6) = E1; K(7) = E2;

  options = odeset('RelTol',RelTol); 
  [t, SOL] = ode45(@(t,V) FNode(t,V,K), tSpan, V0, options);

% Membrane potential  [mV]
  Vm = SOL(:,1);
  Vinf = R.*Iext + E1;
% Resistor current [ohms]
  IR = (Vm-E1)./R;
  IC = Iext - IR;

% GRAPHICS  ===========================================================

figure(1)
   pos = [0.05 0.05 0.29 0.69];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
subplot(3,1,1)   
xP = tm; yP = Vm.*1e3;   
   plot(xP,yP,'b','linewidth',2)
   grid on
   ylabel('V_M  [ mV ]')
   xlabel('t  [ ms ]')
   set(gca,'fontsize',12)
 
subplot(3,1,2) 
   xP = tm; yP = E1.*ones(N,1).*1e3;
   plot(xP,yP,'r','linewidth',1)
   hold on
   xP = tm; yP = Vinf.*1e3;   
   plot(xP,yP,'b','linewidth',2)
     
   grid on
   ylabel('V_{inf}  [ mV ]')
   xlabel('t  [ ms ]')
   set(gca,'fontsize',12)   

subplot(3,1,3)   
   xP = tm; yP = -Iext.*1e9;   
   plot(xP,yP,'b','linewidth',2)
   hold on
   xP = tm; yP = IR.*1e9;   
   plot(xP,yP,'r','linewidth',2)
   xP = tm; yP = IC.*1e9;   
   plot(xP,yP,'m','linewidth',2)
   grid on
   ylabel('I_{ext}  [ nA ]')
   xlabel('t  [ s ]')
   set(gca,'fontsize',12)
      

toc
   
% FUNCTIONS  ==========================================================

function dV = FNode(t,V,K)
  C = K(1);R = K(2);
  tON = K(3); tOFF = K(4); Imax = K(5); E1 = K(6); E2 = K(7);

  Iext = 0;
  if t > tON*1e-3; Iext = Imax; end
  if t > tOFF*1e-3; Iext = 0;     end
 
  dV = -(V-E1-E2)/(R*C) + Iext/C;
 
end