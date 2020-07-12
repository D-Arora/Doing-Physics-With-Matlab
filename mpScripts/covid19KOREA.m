% covid19KOREA.m

% Revised:  2020 05 16

% Mathemtatical model for the spread of the covid19 virus
% Model for SOUTH KOREAL spread of virus
% Data from https://www.worldometers.info/coronavirus/
% Model specified by the variables in the INPUT SECTION of the Script

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/sird19EA.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b

% Populations: Data d / Model
% Susceptible individuals:  Sd  S
% Active infections (currently infected);  Id  I
% Removed: Rm
% Recovered:  Rd  R
% Dead:  Dd  D


clear 
clc
close all

global Imax Rmax Dmax ITmax Smax


% MODEL  ============================================================
% Initialize arrays
  num = 5000;                    % Model: number of time steps 
  S  = zeros(num,1); S(1) = 1;   % Susceptible population
  I = zeros(num,1);              % Active infected population
  Rm = zeros(num,1);             % Removals form infected population  

  
% INPUTS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
% Maximum number of days for model   [200]
  tMax = 200;
% Initial number of infections   [5e-2] 
  I(1) = 5e-2;
% Population scaling factor   [1e4]
  f = 1.05e4;
% S --> I  [0.40] 
   a = 0.40;
% I --> Rm  [0.039]
   b = 0.035;
% Fraction of recoveries compared to deaths   [0.973]   
   r = 0.974;     % recoveries
% Starting Date
  zSTART = datenum([2020 2 26]);
   
% SETUP ===============================================================   
% Fraction of deaths
  d = 1 - r;    
% TIME [days]
  t = linspace(0,tMax,num)';
  dt = t(2) - t(1);  
% Time index
  t_index = 50;  
   
% TIME EVOLUTION OF POPULATIONS  ======================================
  Rm(1) = 200/f;
for n = 1:num-1
  dS = -a*S(n)*I(n)*dt;
  dI = -dS;
  dRm = b*I(n)*dt;
 
  S(n+1) =  S(n)  + dS;
  I(n+1) =  I(n)  + dI - dRm;
  Rm(n+1) = Rm(n) + dRm;
  
  if n == 2220; S(n+1) = 0.2; end
  
end
  
  I = f.*I;
  Rm = f.*Rm;
  D0 = 400;
  k0 = 1.1e-4;
  D = D0.*(1 - exp(-k0.*Rm));
  R = Rm - D;
  Itot = I + Rm;
 
  Imax = max(I);
  Rmax = max(R);
  Dmax = max(D);
  Smax = f*max(S);
  ITmax = max(Itot);  

% Percentages: model
  pI = 100*I(end) / Itot(end);
  pR = 100*R(end) / Itot(end);
  pD = 100*D(end) / Itot(end);
  
  
% DATA  ===============================================================
% Load data from sir function
%  Col 1 Total Cases  Idtot / Col 2 Active Cases Id / Col 3 Deaths Dd
  covid = sir;
  Ndays = length(covid);
  Idtot = covid(:,1);
  Id    = covid(:,2);
  Dd    = covid(:,3);
  Rmd    = Idtot - Id - Dd;    % Removals
  
  Rd = Rmd; 
  Idmax = max(Id);
  Rdmax = max(Rd);
  Ddmax = max(Dd);
  ITdmax = max(Idtot);  
    
% Percentages: Data
  pId = 100*Id(end) / Idtot(end);
  pRd = 100*Rd(end) / Idtot(end);
  pDd = 100*Dd(end) / Idtot(end);
  
% Exponential growth at begiing of virus spread
  th = 4.4; t0 = 0; tE = 0 : 40;     %Ndays;
  k = log(2)/th;
  v0 = Idtot(1);
  v = v0.*exp(k*(tE - t0));   

  
% GRAPHICS  ===========================================================

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.02 0.46 0.80]);
  set(gcf,'color','w');

subplot(3,2,1)   % Total infections
    xP = tE;
    yP = v;
    plot(xP,yP,'r','linewidth',1.2)
   
    hold on
    xP = 1: Ndays;
    yP = Idtot; 
    plot(xP,yP,'b+')
    hold on
    xP = t; yP = Itot;
    plot(xP,yP,'b','linewidth',2)
    
    ylim([0 1.2*ITmax])
    grid on
    xlabel('days elapsed')
    ylabel('I_{tot}','FontName','Times New Roman')
    title('Total Infections I_{tot}','fontweight','normal');
    set(gca,'fontsize',12)  
    set(gca,'Position', [0.0850 0.7093+0.04 0.3347 0.2157]);
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    
  
subplot(3,2,2)  % INFECTIONS
    xP = 1: Ndays;
    yP = Id;
    plot(xP,yP,'b+')   
    hold on
    
%     xP = t_index; yP = Id(t_index);
%     plot(xP,yP,'bo','markerfacecolor','b','markersize',8)
    
    xP = t; yP = I;
    plot(xP,yP,'b','linewidth',2)
   
    grid on
    title('Active Infections  I','fontweight','normal')
    ylabel('I','FontName','Times New Roman') 
    xlabel('days elapsed')
    set(gca,'fontsize',12)
    set(gca,'Position', [0.5703-0.04 0.7093+0.04 0.3347 0.2157]);
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    
    
subplot(3,2,3)   % RECOVERIES
    xP = 1: Ndays;
    yP = Rd; 
    plot(xP,yP,'k+')
    hold on
 
%     xP = t_index; yP = Rd(t_index);
%     plot(xP,yP,'ko','markerfacecolor','k','markersize',8)
%     
    xP = t; yP = R;
    plot(xP,yP,'k','linewidth',2)
    
    grid on
    ylabel('R')
    xlabel('days elapsed')
    title('Recoveries  R','fontweight','normal')
    set(gca,'fontsize',12)
    set(gca,'Position', [0.0800 0.4096+0.04 0.3347 0.2157]);
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    
    
 subplot(3,2,4)    % DEATHS
    xP = 1: Ndays;
    yP = Dd; 
    plot(xP,yP,'r+')
    hold on
 
%     xP = t_index; yP = Dd(t_index);
%     plot(xP,yP,'ro','markerfacecolor','r','markersize',8)
    
    xP = t; yP = D;
    plot(xP,yP,'r','linewidth',2)
   
    grid on
    ylabel('D')
    xlabel('days elapsed')
    title('Deaths  D','fontweight','normal')
    set(gca,'fontsize',12) 
    set(gca,'Position', [0.5703-0.04 0.4096+0.04 0.3347 0.2157]);
    yMax = ylim;
    z = yMax(2);
    plotMonths(z) 
 

 subplot(3,2,5)
    xP = t; yP = S;
    plot(xP,yP,'b','linewidth',2)
    hold on
    
    xlim([0 tMax])
    ylim([0 1.1*Smax/f])
    grid on
    xlabel('days elapsed')
    ylabel('S') 
    title('Susceptible Population S','fontweight','normal')
    set(gca,'fontsize',12)    
    set(gca,'Position', [0.0800 0.1100+0.04 0.3347 0.2157]); 
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
 
subplot(3,2,6)  
   xlim([0 120])
   ylim([0 250])
   hh = 250; dh = -28;
   
   z = zSTART;
   z = datetime(z,'ConvertFrom','datenum');
   z.Format = 'dd-MMM-yyyy';
   txt = cellstr(z) ; 
   Htext = text(0,hh,txt,'fontsize',12);
   set(Htext,'color','k')
   
   hh = hh+2*dh;
   txt = 'MODEL PARAMETERS';
   text(0,hh,txt,'fontsize',12)
   
%    hh = hh+dh;
%    txt = sprintf('I(1) = %3.0e   f = %2.2e', I(1)/f, f);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('a = %2.3f   b = %2.3f   r = %2.3f  d = %2.3f', ...
%        a,b,r,d);
%    text(0,hh,txt,'fontsize',12)
   hh = hh+dh;
   txt = sprintf('I(1) = %3.0e   f = %2.2e   a = %2.3f   b = %2.3f', I(1)/f, f,a ,b);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+2*dh; 
   txt = 'DATA';
   text(0,hh,txt,'fontsize',12)
   
   z = zSTART;
   z = z + Ndays-1;
   z = datetime(z,'ConvertFrom','datenum');
   z.Format = 'dd-MMM-yyyy';
   txt = cellstr(z) ; 
   text(20,hh,txt,'fontsize',12)
    
   hh = hh+dh;
   z1 = Id(end); z2 = Rd(end); z3 = Dd(end); z4 = z1+z2+z3;  
   txt = sprintf('I = %5.0f    R = %5.0f    D = %5.0f    I_{tot} = %5.0f', z1,z2,z3,z4);
   text(6,hh,txt,'fontsize',12)
   
   % Recoveries and Deaths
     hh = hh+dh;
     txt = sprintf('percent: Active = %2.1f  Recoveries = %2.1f    Deaths = %2.1f   ',pId,pRd,pDd);
     text(2,hh,txt,'fontsize',12)  
   
     hh = hh+2*dh; 
     txt = 'MODEL PREDICTIONS';
     text(0,hh,txt,'fontsize',12)
   
     z = zSTART;
     z = z + tMax;
     z = datetime(z,'ConvertFrom','datenum');
     z.Format = 'dd-MMM-yyyy';
     txt = cellstr(z) ; 
     text(78,hh,txt,'fontsize',12)
   
     hh = hh+dh;
     z1 = I(end); z2 = R(end); z3 = D(end); z4 = z1+z2+z3;  
     txt = sprintf('I = %5.0f    R = %5.0f    D = %5.0f    I_{tot} = %5.0f', z1,z2,z3,z4);
     text(6,hh,txt,'fontsize',12)

     hh = hh+dh;
     txt = 'Peak';
     text(10,hh,txt,'fontsize',12)
     z = find(I == Imax,1);
     z = t(z) + zSTART;
     z = datetime(z,'ConvertFrom','datenum');
     z.Format = 'dd-MMM-yyyy';
     txt = cellstr(z) ; 
     text(30,hh,txt,'fontsize',12)
   
     txt = sprintf('I_{peak} = %5.0f',Imax);
     text(80,hh,txt,'fontsize',12)
   
     hh = hh+dh;
     txt = sprintf('percent: Active = %2.1f  Recoveries = %2.1f    Deaths = %2.1f   ',pI,pR,pD);
     text(2,hh,txt,'fontsize',12)   
        
     axis off
     set(gca,'Position', [0.5703-0.08 0.1100+0.04 0.3347 0.2157]);
     

 figure(2)  % ----------------------------------------------------------
   set(gcf,'units','normalized');
   set(gcf,'position',[0.460 0.05 0.2 0.25]);
   set(gcf,'color','w');

   xP = Rd+Dd; yP = Dd;
   plot( xP,yP,'b+')
   hold on
   xP = Rm; yP = D;
   plot( xP,yP,'r','linewidth',2)
   grid on
   ylabel('deaths D')
   xlabel('removals R_m ')
%   title('D = D_0 [1 - exp( - k R_m)]','fontweight','normal')
   txt = sprintf('D_0 = %2.0f', D0);
   text(7000,120,txt,'fontsize',12)
   txt = sprintf('k = %2.2e',k0);
   text(7000,70,txt,'fontsize',12)
   set(gca,'fontsize',12) 
   
     

figure(3)  % ----------------------------------------------------------
   set(gcf,'units','normalized');
   set(gcf,'position',[0.670 0.05 0.2 0.50]);
   set(gcf,'color','w');

subplot(2,1,1)   
   yP = Rd+Dd; xP = Idtot;
   plot( xP,yP,'b+')
   hold on
   yP = Rm; xP = Itot;
   plot( xP,yP,'r','linewidth',2)
   grid on
   xlabel('infections   I_{tot}')
   ylabel('removals   R_m')
   set(gca,'fontsize',12)   

 subplot(2,1,2)   
   yP = Rd+Dd; xP = Id;
   plot( xP,yP,'b+')
   hold on
   yP = Rm; xP = I;
   plot( xP,yP,'r','linewidth',2)
   grid on
   xlabel('active infections   I')
   ylabel('removals   R_m')
   set(gca,'fontsize',12)   
   

figure(9)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.46 0.4 0.2 0.50]);
  set(gcf,'color','w');

  subplot(3,1,1)
    xP = Rd+Dd; yP = Dd;
   plot( xP,yP,'b+')
   hold on
   xP = Rm; yP = D;
   plot( xP,yP,'r','linewidth',2)
   grid on
   ylabel('deaths D')
   xlabel('removals R_m ')
 %  title('D = D_0 [1 - exp( - k R_m)]','fontweight','normal')
   txt = sprintf('D_0 = %2.0f', D0);
   text(7000,120,txt,'fontsize',12)
   txt = sprintf('k = %2.2e',k0);
   text(7000,70,txt,'fontsize',12)
   set(gca,'fontsize',12)
   
 subplot(3,1,2)  
   yP = Rd+Dd; xP = Id;
   plot( xP,yP,'b+')
   hold on
   yP = Rm; xP = I;
   plot( xP,yP,'r','linewidth',2)
   grid on
   xlabel('active infections   I')
   ylabel('removals   R_m')
   set(gca,'fontsize',12)   

 subplot(3,1,3)  
   yP = Rd+Dd; xP = Idtot;
   plot( xP,yP,'b+')
   hold on
   yP = Rm; xP = Itot;
   plot( xP,yP,'r','linewidth',2)
   grid on
   xlabel('infections   I_{tot}')
   ylabel('removals   R_m')
   set(gca,'fontsize',12)   

      
   
% FUNCTIONS ============================================================
  
function covid = sir
c = 1;
T = [1261 1225 12];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1766 1729 13];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2337 2297 16];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [3150 3109 17];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [3736 3685 21];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [4335 4277 28];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [5186 5120 32];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [5621 5498 35];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [6284 6107 42];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [6593 6415 43];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7041 6875 48];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7313 7097 50];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7478 7178 53];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7513 7165 60];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7755 7362 60];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7869 7293 66];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7979 7198 67];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8086 7180 72];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8162 7253 75];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8236 7024 75];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8320 6838 81];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8413 6789 84];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8565 6257 91];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8652 6325 94];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8799 6085 102];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8897 5884 104];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8961 5684 111];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9031 5410 120];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9241 5281 126];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9241 4966 131];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [9332 4665 139];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9478 4523 144];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9583 4398 152];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9661 4275 158];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9786 4216 162];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9887 4155 165];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9976 3979 169];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10062 3867 174];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10156 3654 177];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10237 3591 183];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10284 3500  186];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10331 3445 192];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10384 3408 200];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10423 3246 204];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10450 3125 208];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10480 3026 211];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10512 2930 214];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10537 2873 217];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10564 2808 222];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10591 2750 225];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10613 2627 229];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10635 2576 230];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10653 2484 232];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10661 2385 234];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10674 2324 236];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10683 2233 237];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10694 2179 238];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10702 2051 240];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10708 1967 240];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
   T = [10718 1843 240];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10728 1769 242];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10738 1731 243];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10752 1654 244];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10761 1593 246];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10765 1459 247];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10774 1454 248];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10780 1407 250];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10793 1360 250];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10801 1332 252 252];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10804 1267 254];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10806 1218 255];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10870 1135 256];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10822 1082 256];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10840 1016 256];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10874 1008 256];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10909 1021 256];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10936 1008 258];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10962 1008 259];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10991 969 260];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11018 937 260];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11037 924 262];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11065 898 263];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11065 898 263];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11037 924 262];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11122 723 264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11122 723 264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11142 716 264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11165 705 266];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11190 711 266];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11206 713 267];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11225 681 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11265 701 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11344 735 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11402 770 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11441 774 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11468 793 270];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11503 810 271];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11541 823 272];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11590 850 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11629 857 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11668 889 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11719 915 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11776 951 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11852 989 274];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11902 1015 276];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11902 1017 276];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11947 1017 276];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;



end





function  plotMonths(z)
    
       yMax = z;    
     % Months
       m1 = [2020 02 26];    % day zero
       m2 = [2020 03  1];    % 
       m3 = [2020 04  1];    % 
       m4 = [2020 05  1];    % 
       m5 = [2020 06  1];    % 
       m6 = [2020 07  1];    % 
       m7 = [2020 08  1];    % 
       m8 = [2020 09  1];    % 
   
      tm(1) = datenum(m1);
      tm(2) = datenum(m2) - tm(1);
      tm(3) = datenum(m3) - tm(1);
      tm(4) = datenum(m4) - tm(1);
      tm(5) = datenum(m5) - tm(1);
      tm(6) = datenum(m6) - tm(1);
      tm(7) = datenum(m7) - tm(1);
      tm(8) = datenum(m8) - tm(1);
   
      for n = 2:8
        plot([tm(n), tm(n)], [0, yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
        hold on
      end
      
      ys = 0.9;
      M = 'MAMJJAS';
      for c = 1:6
        text(5+tm(c+1), ys*yMax,M(c))
      end

end
