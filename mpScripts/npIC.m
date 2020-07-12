% mp/mpscripts/npIC.m
 
% ION CHANNELS and GATE VARIABLES
% A finite difference method is used to solve the
%   ordinary differential equation relating to ion channels and gate variables
%   to membrane currents

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/npIC.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  MatlabVisualPhysics@gmail.com
% 191126


clear
close all
clc


% MODEL PARAMETERS ====================================================

% Resting value for membrane potential  [-65e-3  V]
  Vr = -65e-3;
% Membrane potential  [-100e-3   50e-3  V]
  Vmin = -100e-3;   Vmax = 50e-3;
% Reversal potentials [Na+  50e-3    K+ -77e-3   V] 
   VNa = 50e-3;       
   VK = -77e-3;       
% Temperature  [20 degC]
  T = 20;
% Number of grid points  ( N such that dV = Vm - Vr <> -10e-3 ) 
  N = 510;
  
% SETUP ============================================================

% Membrane potential  [V]
  Vm = linspace(Vmin,Vmax,N);
% Initialize arrays
  alphaN = zeros(N,1);
  betaN  = zeros(N,1);
  alphaM = zeros(N,1);
  betaM  = zeros(N,1);
  alphaH = zeros(N,1);
  betaH  = zeros(N,1);
% LinwWidth for plots
  LW = 2;
% FontSize
  FS = 12;
% Scaling Factor  V to mV
  sf = 1e3;
  
  
% RATE CONSTANTS alpha and beta  ===================================  

for c = 1 : N
    [An, Bn, Am, Bm, Ah, Bh] = AlphaBeta(Vm(c), Vr, T); 
    alphaN(c) = An;
    betaN(c)  = Bn;
    alphaM(c) = Am;
    betaM(c)  = Bm;
    alphaH(c) = Ah;
    betaH(c)  = Bh;
end

% Gating variables m,h and n: asymptotic (steady-state) values
%    time constants  [ms]
  nINF = alphaN./(alphaN + betaN);
  tauN = 1./(alphaN + betaN);

  hINF = alphaH./(alphaH + betaH);
  tauH = 1./(alphaH + betaH);
  
  mINF = alphaM./(alphaM + betaM);
  tauM = 1./(alphaM + betaM);
  
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.4 0.7]);
  set(gcf,'color','w');
  
  subplot(3,2,1)
  xP = Vm.*sf; yP = alphaM;
  plot(xP,yP,'r','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\alpha_m','fontsize',18)
  title('sodium','fontweight','normal')
  
  subplot(3,2,3)
  xP = Vm.*sf; yP = betaM;
  plot(xP,yP,'r','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\beta_m','fontsize',18)
  title('sodium','fontweight','normal')
  
 subplot(3,2,2)
  xP = Vm.*sf; yP = alphaN;
  plot(xP,yP,'b','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\alpha_n','fontsize',18)
  title('potasium','fontweight','normal') 
  
  subplot(3,2,4)
  xP = Vm.*sf; yP = betaN;
  plot(xP,yP,'b','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\beta_n','fontsize',18)
  title('potasium','fontweight','normal') 
  
  subplot(3,2,5)
  xP = Vm.*sf; yP = alphaH;
  plot(xP,yP,'r','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\alpha_h','fontsize',18)
  title('sodium','fontweight','normal') 
  
  subplot(3,2,6)
  xP = Vm.*sf; yP = betaH;
  plot(xP,yP,'r','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\beta_h','fontsize',18)
  title('sodium','fontweight','normal')

  % ===================================================================
  figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.43 0.05 0.4 0.7]);
  set(gcf,'color','w');
  
  subplot(3,2,1)
  xP = Vm.*sf; yP = mINF;
  plot(xP,yP,'r','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('m_{INF}','fontsize',18)
  title('sodium','fontweight','normal')
  
  subplot(3,2,2)
  xP = Vm.*sf; yP = tauM;
  plot(xP,yP,'r','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\tau_m  [ ms ]','fontsize',18)
  title('sodium','fontweight','normal')
  
 subplot(3,2,3)
  xP = Vm.*sf; yP = hINF;
  plot(xP,yP,'r','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('h_{inf}','fontsize',18)
  title('sodium','fontweight','normal') 
  
  subplot(3,2,4)
  xP = Vm.*sf; yP = tauH;
  plot(xP,yP,'r','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\tau_h  [ ms ]','fontsize',18)
  title('sodium','fontweight','normal') 
  
  subplot(3,2,5)
  xP = Vm.*sf; yP = nINF;
  plot(xP,yP,'b','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('n_{inf}','fontsize',18)
  title('potasium','fontweight','normal') 
  
  subplot(3,2,6)
  xP = Vm.*sf; yP = tauN;
  plot(xP,yP,'b','linewidth',LW);
  grid on
  box on
  set(gca,'fontsize',FS);
  xlabel('V_m  [ mV ]')
  ylabel('\tau_n  [ ms ]','fontsize',18)
  title('potasium','fontweight','normal') 
  
  
  
% Sodium
% FUNCTIONS =========================================================
  function [An, Bn, Am, Bm, Ah, Bh] = AlphaBeta(V,Vr,T)
   

   V = V*1000;
   Vr = Vr*1000;
   dV = (V - Vr);
   phi = 3^((T-6.3)/10);

   An = phi * (eps + 0.10 - 0.01 .* dV) ./ (eps + exp(1 - 0.1 .* dV) - 1);
   Am = phi * (eps + 2.5 - 0.1 .* dV) ./ (eps + exp(2.5 - 0.1 .* dV) - 1);
   Ah = phi * 0.07 .* exp(-dV ./ 20);

   Bn = phi * 0.125 .* exp(-dV ./ 80);
   Bm = phi * 4 .* exp(-dV/18);
   Bh = phi * 1 ./ (exp(3.0 - 0.1 .* dV) + 1);

  end  