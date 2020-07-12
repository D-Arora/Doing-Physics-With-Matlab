% CN04.m

% RC CIRCUIT FOR STROBE LIGHTING
% Finite Difference Method for RC circuit anlysis

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
 
% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm
% 180216

clear all
close all
clc

% INPUTS ----------------------------------------------------------------
  R = 9e3;
  C = 666e-6;

  VS = 600;
  VB = 300;
  Vtube = 50;
  N = 200;



% CALCULATIONS  ---------------------------------------------------------
  tau = R*C;
  tMax = 5*tau;

  t = linspace(0,tMax, N);
  dt = t(2)-t(1);

% Initialize arrays
  vS = zeros(N,1);
  vR = zeros(N,1);
  vC = zeros(N,1);
  iS = zeros(N,1);

% Specify source emf
  N1 = round(N/20);
  vS(N1:N) = VS;

% Finite difference algorithm -------------------------------------------
   k = dt/C;
% Time Steps #2 to #N
  for c = 2 : N
    vC(c)  = vC(c-1) + k*iS(c-1);
    vR(c)  = vS(c) - vC(c);
    iS(c)  = vR(c)/R;
    if vC(c) > VB; vC(c) = Vtube; end
  end

% Calculate flash rate ------------------------------------------------
  [iS_Peaks, t_Peaks] = findpeaks(real(vC));
  nPeaks = length(t_Peaks)-1;
  T_Peaks = (t(t_Peaks(end)) - t(t_Peaks(end-nPeaks)))/(nPeaks);
  f_Peaks = 1 / T_Peaks;  
  
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
  tm1 = 'T_{peaks} = ';
  tm2 = num2str(T_Peaks,'%3.2f');
  tm3 = '  s     ';
  tm4 = 'f_{peaks} = ';
  tm5 = num2str(f_Peaks,'%3.2f');
  tm6 = '  Hz';
  tm = [tm1 tm2 tm3 tm4 tm5 tm6];
  title(tm);
  set(gca,'yLim',[0 650]);
  set(gca,'fontsize',12);


subplot(2,1,2)
  plot(t,vC,'r','linewidth',2)
  ylabel('v_C  [ V ]');
  xlabel('t [ s ]');
  grid on
  set(gca,'fontsize',12);
  
  