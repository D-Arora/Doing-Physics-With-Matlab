% mp/mpscripts/npHHA.m
 
% HODGKIN-HUXLEY MODEL FOR MEMBRANE CURRENT
% A finite difference method is used to solve the
%   ordinary differential equation relating membrane potential
%   to membrane currents

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/??
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  MatlabVisualPhysics@gmail.com
% 191113



close all
clear 
clc

global Vr

tic
% INPUT: INJECTED (EXTERNAL)) CURRENT =================================
% flagJ = 1 2 3 4 5 6

  flagJ = 1;

% Variables set in switch / case statements

 
 
switch flagJ
     
    case 1  % Single current pulse
        
      % Amplitude of pulse  [100e-6  A.cm-2] 
        J0 = 100e-6;
      % Simulation time [5.0e-3  s] 
        tMax = 10e-3;
      % Time stimulus applied [0.5e-3  s] 
        tStart = 0.5e-3;    
      % Pulse duration: ON time  [0.1e-3 s]
        tON = 0.1e-3;
      % Number of grid points  [8001] 
        num = 8001;  
      
      % external current density [A.cm^-2]
        Jext = zeros(num,1);        
        t = linspace(0,tMax,num);
        dt = t(2) - t(1);  
        tOFF = tON + tStart;
        num1 = find(t > tStart, 1 );    % index for stimulus ON
        num2 = find(t > tOFF, 1 );   % index for stimulus OFF
        Jext(num1:num2) = J0;        % external stimulus current
      
      % Stimulus: electric charge per unit area
      Qext = J0 * tON;
      disp('Stimulus: Electric charge injected')
      fprintf(' Qext = %2.2e  C/cm^-2 \n',Qext)
      disp('  ')
      
      
    case 2    % Constant current injection step function
       % Amplitude of step function  [100e-6  A.cm-2] 
        J0 = 10e-6;
      % Simulation time [40e-3  s] 
        tMax = 40e-3;
      % Time stimulus applied [5e-3  s] 
        tStart = 5e-3;    
      % Number of grid points  [8001] 
        num = 8001;     
        
        Jext = zeros(num,1);       % external current density [A.cm^-2] 
        t = linspace(0,tMax,num);
        dt = t(2) - t(1);  
        num1 = find(t > tStart, 1 );   % index for stimulus ON
        Jext(num1:num) = J0;     % external stimulus current 
    
        
    case 3    % Double current pulse stimulus
      % Amplitude of pulses  [100e-6  A.cm-2] 
        J0 = 100e-6;
      % Simulation time [2.5e-3  s] 
        tMax = 10e-3;
      % Pulse duration: ON time  [0.1e-3 s]
        tON = 0.1e-3;
      % time pulse #1 ON  [1.0e-3 s] / time pulse #2 ON  [5.0e-3 s]
        t1 = 1.0e-3;
        t2 = 3.0e-3; 
      % Number of grid points  [8001] 
        num = 8001;       
        
        Jext = zeros(num,1);       % external current density [A.cm^-2] 
        t = linspace(0,tMax,num);
        dt = t(2) - t(1);  
        num1 = find(t > t1,1);   % index for stimulus #1 ON
        num2 = find(t > t1+tON, 1 );   % index for stimulus #1 OFF
        Jext(num1:num2) = J0;    % external stimulus current
        num1 = find(t > t2,1);   % index for stimulus #2  ON
        num2 = find(t > t2+tON, 1 );   % index for stimulus #2 OFF
        Jext(num1:num2) = J0;    % external stimulus current
      
       
    case 4  % Series of current pulses
      % Amplitude of pulses [100e-6  A.cm-2] 
        J0 = 100e-6;
      % Simulation time [2.5e-3  s] 
        tMax = 25e-3; 
      % Stimulus applied  [1e-3  s]
        tSTART = 1e-3; 
      % Pulse ON / OFF times  [1e3  2e3]
        tON = 1e-3;
        tOFF = 0.5e-3;
      % Number of grid points  [8001] 
        num = 8001;   
        
        Jext = zeros(num,1);       % external current density [A.cm^-2] 
        t = linspace(0,tMax,num);
        dt = t(2) - t(1);     
        
      % Number of pulses
        tPulse = tON+tOFF; 
        nPulse = (t(end) - tSTART) / (tON + tOFF);
      
     for cn = 1:floor(nPulse)
       t1 = tSTART + (cn-1)* tPulse;
       t2 = t1 + tON;
       num1 = find(t > t1,1);
       num2 = find(t > t2,1);
       Jext(num1:num2) = J0;    % external stimulus current
     end
         
      
    case 5   % Sinusoidal current stimulus
      % Amplitude of step function  [10e-6  A.cm-2] 
        J0 = 10e-6;
      % Simulation time [50e-3  s] 
        tMax = 50e-3; 
      % Frequency of stimulus  [200 Hz]
        pf = 400; 
      % Number of grid points  [1001] 
        num = 1001;     
               
        pT = 1 / pf; 
        
        t = linspace(0,tMax,num);
        dt = t(2) - t(1);     
        
        Jext = J0 .* sin(2*pi*t/pT);
   
          
    case 6  % Noise
       % Amplitude of step function  [7e-6  A.cm-2] 
        J0 = 7e-6;
      % Simulation time [2.5e-3  s] 
        tMax = 50e-3; 
      % Stimulus applied    [1e-3  s]  
        tSTART = 1e-3; 
      % Pulse ON / OFF times
        tON = 1e-3;
        tOFF = 2e-3;
      % Number of grid points  [8001] 
        num = 8001;   
        
        Jext = zeros(num,1);       % external current density [A.cm^-2] 
        t = linspace(0,tMax,num);
        dt = t(2) - t(1);     
        
      % Number of pulses
        tPulse = tON+tOFF; 
        nPulse = (t(end) - tSTART) / (tON + tOFF);
      
        rng('shuffle');
     for cn = 1:floor(nPulse)
       t1 = tSTART + (cn-1)* tPulse;
       t2 = t1 + tON;
       num1 = find(t > t1,1);
       num2 = find(t > t2,1);
       Jext(num1:num2) = J0 + 3*(J0) .* (2.*rand(1,1)-1);    
     end
     
end


% FIXED PARAMETERS ====================================================

   sf = 1e3;          % conversion of V to mV
   VR = -65e-3;       % resting voltage [V]
   Vr = VR*1e3;       % resting voltage [mV]
   VNa = 50e-3;       % reversal voltage for Na+ [V]
   VK = -77e-3;       % reversal voltage for K+  [V]
   Cm = 1e-6;         % membrane capacitance/area  [F.cm^-2]

   gKmax = 36e-3;     % K+ conductance [S.cm^-2]
   gNamax = 120e-3;   % Na+ conductance [S.cm.-2)]
   gLmax = 0.3e-3;    % max leakage conductance [S.cm-2]
   T = 20;          % temperature [20 deg C]


% SETUP ===============================================================

  JNa  = zeros(num,1);       % Na+ current density (A.cm^-2)
  JK   = zeros(num,1);       % K+  current density (A.cm^-2)
  JL   = zeros(num,1);       % leakage current density (A.cm^-2)
  Jm   = zeros(num,1);       % membrane current (A.cm^-2)
  V    = zeros(num,1);       % membrane potential (V)
  gNa  = zeros(num,1);       % Na+ conductance
  gK   = zeros(num,1);       % K+ conductance
  gL   = ones(num,1);        % gL conductance
  n    = zeros(num,1);       % K+ gate parameter
  m    = zeros(num,1);       % Na+ gate parameter
  h    = zeros(num,1);       % Na+ gate parameter


% Initial Values  -----------------------------------------------------
  V(1) = VR;                   % initial value for membrane potential
  [An, Bn, Am, Bm, Ah, Bh] = AlphaBeta(V(1), T);

  n(1) = An / (An + Bn);
  m(1) = Am / (Am + Bm);
  h(1) = Ah / (Ah + Bh);

  gK(1)  = gKmax * n(1)^4;
  gNa(1) = gNamax * m(1)^3 * h(1);
  gL = gLmax .* gL;

  JK(1)  = gK(1)  * (V(1) - VK);
  JNa(1) = gNa(1) * (V(1) - VNa);

  Vadj = (Jext(1) - JK(1) - JNa(1))/gL(1);
  JL(1)  = gL(1) * (V(1) - VR + Vadj);

  Jm(1)  = JNa(1) + JK(1) + JL(1);

  V(1) = VR + (dt/Cm) * (-JK(1) - JNa(1) - JL(1) + Jext(1));


% Solve ODE -----------------------------------------------------------

for cc = 1 : num-1
    
 [An, Bn, Am, Bm, Ah, Bh] = AlphaBeta(V(cc), T);

  An = sf * An;   Am = sf * Am;   Ah = sf * Ah;  
  Bn = sf * Bn;   Bm = sf * Bm;   Bh = sf * Bh; 

  n(cc+1) = n(cc) + dt * (An *(1-n(cc)) - Bn * n(cc)); 
  m(cc+1) = m(cc) + dt * (Am *(1-m(cc)) - Bm * m(cc)); 
  h(cc+1) = h(cc) + dt * (Ah *(1-h(cc)) - Bh * h(cc)); 

  gK(cc+1) = n(cc+1)^4 * gKmax;
  gNa(cc+1) = m(cc+1)^3 * h(cc+1) * gNamax;

  JK(cc+1)  = gK(cc+1)  * (V(cc) - VK);
  JNa(cc+1) = gNa(cc+1) * (V(cc) - VNa);
  JL(cc+1)  = gL(cc+1) * (V(cc) - VR + Vadj);
  Jm(cc+1)  = JNa(cc+1) + JK(cc+1) + JL(cc+1);

  V(cc+1) = V(cc) + (dt/Cm) * (-JK(cc+1) - JNa(cc+1) - JL(cc+1) + Jext(cc+1));

end

 Vdot = (Jext - Jm)./Cm;

% GRAPHICS ============================================================

% Width & height  of Figure Window
  w = 0.3; d = 0.30;
% Position of Figure window (x,y)
  x1 = 0.02; y1 = 0.45;
  x2 = 0.35; y2 = 0.05;
  x3 = 0.67;
  
% Current 
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[x1 y1 w 1.3*d]);
  set(gcf,'color','w');
  LW = 2;

  sf6 = 1e6;
  
  subplot('position',[0.15 0.4 0.8 0.6])
    x = t.*sf; y = JNa.*sf6;
    plot(x,y,'r','linewidth',LW);   %  Current - Na+
  hold on
  y = JK.*sf6;
    plot(x,y,'b','linewidth',LW);   %  Current - K+

  y = JL.*sf6;
    plot(x,y,'c','linewidth',LW);   %  Current - leakage

  y = Jm.*sf6;
    plot(x,y,'k','linewidth',LW);   %  Current - K+

  Hlegend = legend(' J_{Na}',' J_K', 'J_L',' J_m');
    set(Hlegend,'orientation','horizontal','location','northoutside','box','off', ...
        'fontsize',12)
  %xlabel('time  t   [ ms ]');
   ylabel('J  [ \muA.cm ^{-2} ]');
  set(gca,'fontsize',12)
  grid on
  box on
  
 subplot('position',[0.15 0.11 0.8 0.2])
   x = t.*sf;   y = Jext.*sf6;
    plot(x,y,'m','linewidth',LW);       %  Current - ext
  
    xlabel('time  t   [ ms ]'); ylabel('J_{ext}  [ \muA.cm ^{-2} ]');
    set(gca,'fontsize',12)
    grid on
    box on 
    
    
%  Membrane Potential  
figure(2)     
  set(gcf,'units','normalized');
  set(gcf,'position',[x2 y1 w d]);
  set(gcf,'color','w');
  LW = 2;
    
  x = t.*sf;   y = V.*sf;
  plot(x,y,'b','linewidth',LW);   % membrane voltage

  xlabel('time  t   [ ms ]'); ylabel('membrane voltage   V_m  [ mV ]');
  set(gca,'fontsize',12)
  grid on
  box on
  hold on

 if flagJ == 2
   tm = sprintf('J_0  =  %2.0f  \\muA.cm^{-2}  ',J0.*1e6);
   title(tm)
 end
  
%  Conductance
figure(3)     %   
  set(gcf,'units','normalized');
  set(gcf,'position',[x3 y1 w d]);
  set(gcf,'color','w');
  LW = 2;

  x = t.*sf;   y = gNa.*1e3;
   plot(x,y,'r','linewidth',LW);   % conductance  Na+
  hold on
  y = gK.*1e3;
   plot(x,y,'b','linewidth',LW);   % conductance  K+
      
  xlabel('time  t   [ ms ]'); ylabel('conductance  g  [ mS.cm^{-2} ]');
  grid on
  Hlegend = legend(' Na^+ ',' K^+');
    set(Hlegend,'orientation','horizontal','location','northoutside','box','off', ...
        'fontsize',12)
  set(gca,'fontsize',12)
  

%  I-V curve (phase portrait plot)  
figure(4)     
set(gcf,'units','normalized');
  set(gcf,'position',[x2 y2 w d]);
  set(gcf,'color','w');
  LW = 2;
  
 x = V.*sf;   y = Vdot;%Jm.*sf;
  plot(x,y,'b','linewidth',LW); 
  hold on
 
 x = V(1).*sf;   y = Vdot(1); %Jm(1).*sf;
  Hplot = plot(x,y,'go'); 
  set(Hplot,'markersize',8,'markerfacecolor','g');
  hold on 
  
  x = V(end).*sf;   y = Vdot(end); %Jm(end).*sf;
  Hplot = plot(x,y,'ro'); 
  set(Hplot,'markersize',8,'markerfacecolor','r');
  hold on 
  
  xlabel('V_m  [ mV ]');  ylabel('dV_m / dt  [mV/ms]'); %ylabel('J_m  [mA.cm_{-2} ]');
  set(gca,'fontsize',12)
  grid on
  box on  
  
%  Gate Variables n m h 
figure(5)     
  set(gcf,'units','normalized');
  set(gcf,'position',[x3 y2 w d]);
  set(gcf,'color','w');
  LW = 2;

  x = t.*sf;   y = n;
    plot(x,y,'b','linewidth',LW);      
  hold on
  
  y = m;
    plot(x,y,'r','linewidth',LW);   

  y = h;
    plot(x,y,'m','linewidth',LW);   

  Hlegend = legend(' n',' m',' h');
    set(Hlegend,'orientation','horizontal','location','northoutside','box','off', ...
        'fontsize',12)
  xlabel('time  t   [ ms ]'); ylabel('n  m   h');
  set(gca,'fontsize',12)
  grid on
  box on  
  
  
toc
 

% FUNCTIONS =========================================================
  function [An, Bn, Am, Bm, Ah, Bh] = AlphaBeta(V, T)
   global Vr

   V = V*1000;
   dV = (V - Vr);
   phi = 3^((T-6.3)/10);

   An = phi * (eps + 0.10 - 0.01 .* dV) ./ (eps + exp(1 - 0.1 .* dV) - 1);
   Am = phi * (eps + 2.5 - 0.1 .* dV) ./ (eps + exp(2.5 - 0.1 .* dV) - 1);
   Ah = phi * 0.07 .* exp(-dV ./ 20);

   Bn = phi * 0.125 .* exp(-dV ./ 80);
   Bm = phi * 4 .* exp(-dV/18);
   Bh = phi * 1 ./ (exp(3.0 - 0.1 .* dV) + 1);

  end
  
