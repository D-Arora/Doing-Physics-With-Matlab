% CN05.m

% SIMPLE RC CIRCUIT: SQUARE WAVE EXCITATION
% Finite Difference Method for RC circuit anlysis

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
 
% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 180216

clear all
close all
clc

% INPUTS (S.I. Units) ---------------------------------------------------
% resistor [1e4]
  R = 1e4;
  % Capacitor [100e-6]
  C = 10e-6;

% Square wave values VH (high) [10]  VL (low)  [-10]  
  VH = 10;
  VL = -10;
% Number of calculations (check dt << tau)  [5000]   
  N = 5000;
% Max time interval  [100]
  tMax = 100;
% Period of square wave function  [20]  
  T = 20;
  
% CALCULATIONS  ---------------------------------------------------------
  tau = R*C;
  t = linspace(0,tMax, N);
  dt = t(2 )- t(1);

% Initialize arrays
  vS = zeros(N,1);
  vR = zeros(N,1);
  vC = zeros(N,1);
  iS = zeros(N,1);

  % Specify source emf 
  for c = 1 : N
    y = sin(2*pi*t(c)/T);
    if y >= 0
          vS(c) = VH;
    elseif y < 0
          vS(c) = VL;
    end
  end


% Finite difference algorithm -------------------------------------------
   k = dt/C;
% Time Steps #2 to #N
  for c = 2 : N
    vC(c)  = vC(c-1) + k*iS(c-1);
    vR(c)  = vS(c) - vC(c);
    iS(c)  = vR(c)/R;
  end

  
% GRAPHICS -------------------------------------------------------------
  
figure(1)
   pos = [0.07 0.05 0.28 0.47];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
subplot(2,1,1)
   plot(t,vS,'b','linewidth',2)
  grid on
  ylabel('v_S  [ V ]');
  tm1 = '   R = ';
  tm2 = num2str(R,'%3.2e');
  tm3 = ' \Omega  ';
  tm4 = 'C = ';
  tm5 = num2str(C,'%3.2e');
  tm6 = ' F   ';
  tm7 = '\tau = ';
  tm8 = num2str(tau,'%3.2e');
  tm9 = ' s   ';
  tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7 tm8 tm9];
  title(tm);
  set(gca,'yLim',[-12 12]);
  set(gca,'fontsize',12);


subplot(2,1,2)
  plot(t,vC,'r','linewidth',2)
  ylabel('v_C  [ V ]');
  xlabel('t [ s ]');
  grid on
  set(gca,'yLim',[-12 12]);
  set(gca,'fontsize',12);
  
  