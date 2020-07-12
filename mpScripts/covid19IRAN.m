% covid19ITALY.m

% Revised:  2020 05 19

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
  R = zeros(num,1);  
  D = zeros(num,1); 
  
% INPUTS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
% Maximum number of days for model   [200]
  tMax = 200;
% Initial number of infections   [5e-2] 
  I(1) = 1.3e-3;
% Population scaling factor   [1e4]
  f = 10e4;
% S --> I  [0.40] 
   a = 0.25;
% I --> Rm  [0.039]
   b = 0.075;

% Starting Date
  zSTART = datenum([2020 2 26]);
   
% SETUP ===============================================================   

% TIME [days]
  t = linspace(0,tMax,num)';
  dt = t(2) - t(1);  
% Time index
  t_index = 50;  
   
% TIME EVOLUTION OF POPULATIONS  ======================================
  Rm(1) = 44/f;
for n = 1:num-1
  dS = -a*S(n)*I(n)*dt;
  dI = -dS;
  dRm = b*I(n)*dt;
 
  S(n+1) =  S(n)  + dS;
  
   if n == 1700
      a = 0.105;
      b = 0.075;
      S(n+1) = 1.2;
   end
     
   if n == 2200; S(n+1) = 0.9; end
   if n == 2300; S(n+1) = 1.0; end
   
  I(n+1) =  I(n)  + dI - dRm;
  Rm(n+1) = Rm(n) + dRm;
  
end
  
% Scale populations
  I = f.*I;
  Rm = f.*Rm;
  k0 = 1.2e-5;
  D0 = 1e4;
  D = D0.*(1 - exp(-k0.*Rm));
  
  k1 = 2e-5;
  D(1700:end) = D(1700) + 3000.*(1 - exp(-k1.*(Rm(1700:end)-Rm(1700))));
  
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
  th = 5; t0 = -10; tE = 0 : 40;     %Ndays;
  k = log(2)/th;
  v0 = Idtot(1);
  v = v0.*exp(k*(tE - t0));   

  
% GRAPHICS  ===========================================================

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.02 0.435 0.80]);
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
    
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    
  %  ylim([0 1.2*ITmax])
    grid on
    xlabel('days elapsed')
    ylabel('I_{tot}','FontName','Times New Roman')
    title('Total Infections I_{tot}','fontweight','normal');
    set(gca,'fontsize',12)  
    set(gca,'Position', [0.0850 0.7093+0.04 0.3347 0.2157]);
  
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
    yMax = get(gca,'ylim');
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
    ylim([0 1.25])
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
  % title('D = D_0 [1 - exp( - k R_m)]','fontweight','normal')
  % txt = sprintf('D_0 = %2.0f', D0);
  % text(8e4,2.4e3,txt,'fontsize',12)
  % txt = sprintf('k = %2.1e',k0);
  % text(8e4,1.5e3,txt,'fontsize',12)
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
   text(1.7e5,6000,txt,'fontsize',12)
   txt = sprintf('k = %2.2e',k0);
   text(1.7e5,3000,txt,'fontsize',12)
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

T = [139 95 19];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [245 165 26];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [388 281 34];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [593 427 43];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [978 749 54];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [1501 1144 66];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2336 1824 77];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [2922 2278 92];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [3513 2666 108];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [4747 3710 124];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [5823 4009 145];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [6566 4238 194];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [7161 4530 237];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8042 5020 291];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [9000 5687 354];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10075 6370 429];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11364 7321 514];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12729 7779 611];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [13938 8624 724];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14991 9142 853];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [16169 9792 988];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [17361 10516 1135];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [18407 11144 1284];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [19644 11466 1433];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [20610 11419 1556];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [21638 12040 1685];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23049 12861 1812];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [24811 13964 1934];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27017 15315 2077];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [29406 16715 2234];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [32332 18821 2378];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [35408 21212 2517];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [38309 23278 2640];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [41495 24827 2757];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [44605 27051 2898];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [47593 29084 3036];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [50468 30597 3160];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [53183 31954 3294];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [55743 32555 3452];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [58226 32612 3603];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [60500 32525 3739];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [62589 31678 3872];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [64586 30781 3993];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [66220 29801 4110];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [68192 28495 4232];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [70029 23725 4357 ];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [71686 23318 4474];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [73303 22735 4585];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [74877 22065 4683];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [76339 21679 4777];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [77995 20897 4869];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [79494 20472 4958];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [80868 19850 5031];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [82211 20070 5118];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [83505 19023 5209];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [84802 18540 5297];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [85996 17492 5391];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [87026 16702 5481];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [88194 16021 5574];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [89328 15485 5650];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [90481 15114 5710];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [91472 14733 5806];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [92584 14268 5877];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [93657 13909 5957];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [94640 13509 6028];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [95646 13237 6091];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [96448 12942 6156];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [97424 12799 6203];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [98647 12991 6277];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [99970 13155 6340];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [101650 13645 6418];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [103135 13905 6456];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [104691 14313 6541];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [106220 14567 6589];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [107603 14820 6640];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [109286 15179 6685];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [110767 15677 6733];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [112725 16514 6783];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [114533 17140 6854];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [116635 17897 6902];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [118392 18308 6937];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [120198 18746 6988];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [122492 19774 7057];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [124603 20311 7119];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [126949 20958 7183];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [129341 21528 7349];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [131652 22076 7300];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [133521 22090 7359];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [135701 22483 7417];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [137724 22560 7451];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [139511 22566 7508];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [141591 22851 7564];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [143849 23234 7627];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [146668 24060 7677];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [148950 24389 7737];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [151466 24821 7797];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [154445 25563 7878];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [157562 26543 7942 ];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [160696 27478 8012];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [164270 28714 8071];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [167156 29281 8134];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [169425 29178 8209];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [171789 29159 8281];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [173832 29121 8351];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [175927 29045 8425];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [177938 28842 8506];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [180156 28909 8584];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
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
