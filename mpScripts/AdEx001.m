% AdEx001.m


clear
close all
clc


% INPUTS  ================================================================

% V and Vm  membrane potential  [mV}
% w         adaptation current (variable)  [pA]

% membrane capacitance  [pF]
  C = 200;
% synapatic conductance  [nS]  
  gL = 12;
% resting potential (effective)  [mV]  
  EL = -70;
% spike threshold potential  [mV]  
  VT = -50;
% slope factor  [mV]  
  DT = 2;
% level of subthreshold adaptation  [nS]  
  a = 2;
% adaptation time constant  [ms]  
  tw = 300;
% adaptation current increment  [pA]  
  b = 70;
% reset membrane potential  [mV]   
  Vr = -58;
% current step  [pA]  
  I0 = 500;
% threshold potential for action potential  
  Vthres = 0;
% time steps
  N = 5001;
% Max simulation time  
  tMax = 500;



% SETUP  =================================================================

  V = zeros(N,1);
  Vm = V;
  w = V;
  I = zeros(N,1);
  
% Initial values
  V(1) = -58; Vm(1) = V(1);
  w(1) = 0;

  I(500:N) = I0;

% membrane time constat
  tm = C/gL;

% simulation time  
  t = linspace(0,tMax,N);
  dt = t(2) - t(1);
  


% SOLVE ODE  =============================================================

for c = 1 :  N-1

 V(c+1) = V(c) + (dt/C)*( -gL*(V(c)-EL)  + gL*DT*exp( (V(c)-VT)/DT ) +I(c) - w(c) );
 
 w(c+1) = w(c) + (dt/tw)*( a*(V(c+1) - EL) - w(c) ); 
 
 Vm(c+1) = V(c+1);     
 
 if V(c+1) > Vthres
   V(c+1) = Vr; Vm(c+1) = 30;    % spike 
   w(c+1) = w(c+1) + b;
 end
 
end


% NULLCLINES  ===========================================================
  Vn = linspace(-70,-42,501);
 
  w1 =  a.*(Vn - EL);

  w2 = -gL.*(Vn - EL) + (gL*DT).*exp( (Vn-VT)./DT) + I0 ;


% GRAPHICS  ==============================================================
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.02 0.02 0.45 0.65]);
   set(gcf,'color','w');
   FS = 14;
   
   xP = t;
      
   subplot(3,1,1)
   yP = I;
   plot(xP,yP,'b','linewidth',2)
   box on
   grid on
   ylabel('I [pA]')
   set(gca,'fontsize',FS)
   
   subplot(3,1,2)
   yP = Vm;
   plot(xP,yP,'b','linewidth',2)
   box on
   grid on
   ylabel('V_m  [mV]')
   xlabel('t  [ms]')
   set(gca,'fontsize',FS)
   
   
   subplot(3,1,3)
   xP = Vm; yP = w;
   plot(xP,yP,'b','linewidth',2)
   hold on
   hP = plot(xP(1),yP(1),'go');
   set(hP,'markerfacecolor','g')
   box on
   grid on
   
   xP = Vn; yP = w1;
   plot(xP,yP,'m','linewidth',2)
   xP = Vn; yP = w2;
   plot(xP,yP,'r','linewidth',2)
   
   
   ylabel('w  [pA]')
   xlabel('V_m  [mV]')
   set(gca,'fontsize',FS)
   

   
  
   
   