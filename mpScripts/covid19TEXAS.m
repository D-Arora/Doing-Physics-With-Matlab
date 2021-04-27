% covid19TEXAS.m

% Revised:  2020 05 29

% Mathemtatical model for the spread of the covid19 virus
% Model for TEXAS spread of virus
% Data from https://www.worldometers.info/coronavirus/
% Model specified by the variables in the INPUT SECTION of the Script

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2020a

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
  I(1) = 0.0050;
% Population scaling factor   [1e4]
  f = 9e4;
% S --> I  [0.40] 
   a = 0.13;
% I --> Rm  [0.039]
   b = 0.041;
% Fraction of recoveries compared to deaths   [0.973]   
%   r = 0.96;     % recoveries
% Starting Date
  zSTART = datenum([2020 3 12]);
   
% SETUP ===============================================================   
% Fraction of deaths
%  d = 1 - r;    
% TIME [days]
  t = linspace(0,tMax,num)';
  dt = t(2) - t(1);  
% Time index
  t_index = 50;  
   
% TIME EVOLUTION OF POPULATIONS  ======================================
Rm(1) = 28/f;

for n = 1:num-1
  dS = -a*S(n)*I(n)*dt;
  dI = -dS;
  dRm = b*I(n)*dt;
  S(n+1) =  S(n)  + dS;
  
%  if n == 1900; S(n+1) = 0.4; end
%   if n == 2000; S(n+1) = 0.5; end
   if n == 2200; S(n+1) = 0.7; end
   if n == 2400; S(n+1) = 0.8; end
   if n == 2500; S(n+1) = 1; end
   if n == 2600; S(n+1) = 1; end
   if n == 2700; S(n+1) = 1; end
   if n == 2800; S(n+1) = 0.8; end
   if n == 2900; S(n+1) = 0.8; end
   if n == 3000; S(n+1) = 0.8; end
   if n == 3100; S(n+1) = 0.8; end
   
  I(n+1) =  I(n)  + dI - dRm;
  Rm(n+1) = Rm(n) + dRm;
end

% Populations
  I = f.*I;
  Rm = f.*Rm;
  k0 = 2.6e-5;
  D0 = 2.9e3;
  D = D0.*(1 - exp(-k0.*Rm));
  
  k1 = 0.7e-5; D1 = 2000; nS = 2750;
  D(nS:end) = +D(nS) + D1.*(1 - exp(-k1.*(Rm(nS:end)-Rm(nS))));
  
%   D0 = 0.06;
%   D  = D0.*Rm;
%   
  R = Rm - D;
  Itot = I + Rm;
  
  
  
  Imax = max(I);
  Rmax = max(R);
  Dmax = max(D);
  Smax = max(S);
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
  
% Exponential growth at begining of virus spread
  th = 5; t0 = -10; tE = 0 : 60;     %Ndays;
  k = log(2)/th;
  v0 = Idtot(1);
  v = v0.*exp(k*(tE - t0));   

  
% GRAPHICS  ===========================================================

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.1 0.71 0.6])
  set(gcf,'color','w');

subplot(2,3,1)   % Total infections
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
%   set(gca,'Position', [0.0850 0.7093+0.04 0.3347 0.2157]);
    
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    
subplot(2,3,2)  % INFECTIONS
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
 %   set(gca,'Position', [0.5703-0.04 0.7093+0.04 0.3347 0.2157]);
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    
subplot(2,3,3)   % RECOVERIES
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
%    set(gca,'Position', [0.0800 0.4096+0.04 0.3347 0.2157]);
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    
 
 subplot(2,3,4)    % DEATHS
    xP = 1: Ndays;
    yP = Dd; 
    plot(xP,yP,'r+')
    hold on
 
%     xP = t_index; yP = Dd(t_index);
%     plot(xP,yP,'ro','markerfacecolor','r','markersize',8)
    
    xP = t; yP = D;
    plot(xP,yP,'r','linewidth',2)
       
    ylabel('D')
    xlabel('days elapsed')
    title('Deaths  D','fontweight','normal')
    set(gca,'fontsize',12) 
  %  set(gca,'Position', [0.5703-0.04 0.4096+0.04 0.3347 0.2157]);
    
   yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    grid on
    

 subplot(2,3,5)
    xP = t; yP = S;
    plot(xP,yP,'b','linewidth',2)
    hold on
        
    xlim([0 tMax])
  %  ylim([0 1.1])
    grid on
    xlabel('days elapsed')
    ylabel('S') 
    title('Susceptible Population S','fontweight','normal')
    set(gca,'fontsize',12)    
 %   set(gca,'Position', [0.0800 0.1100+0.04 0.3347 0.2157]); 
    
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
    
subplot(2,3,6)  
   xlim([0 120])
   ylim([0 250])
   hh = 250; dh = -28;
   
   txt = 'Start date' ; 
   text(0,hh,txt,'fontsize',12)
   z = zSTART;
   z = datetime(z,'ConvertFrom','datenum');
   z.Format = 'dd-MMM-yyyy';
   txt = cellstr(z) ; 
   Htext = text(35,hh,txt,'fontsize',12);
   set(Htext,'color','k')
   
   hh = hh+dh;
   txt = 'MODEL PARAMETERS';
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('I(1) = %3.0e   f = %2.2e   a = %2.3f   b = %2.3f', I(1)/f,f,a,b);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh; 
   txt = 'DATA';
   text(0,hh,txt,'fontsize',12)
   
   z = zSTART;
   z = z + Ndays;
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
  %   set(gca,'Position', [0.5703-0.08 0.1100+0.04 0.3347 0.2157]);
     

     
%   figure(2)  % ----------------------------------------------------------
%    set(gcf,'units','normalized');
%    set(gcf,'position',[0.460 0.05 0.2 0.25]);
%    set(gcf,'color','w');
% 
%    xP = Rd+Dd; yP = Dd;
%    plot( xP,yP,'b+')
%    hold on
%    xP = Rm; yP = D;
%    plot( xP,yP,'r','linewidth',2)
%    grid on
%    ylabel('deaths D')
%    xlabel('removals R_m ')
%  %  title('D = D_0 [1 - exp( - k R_m)]','fontweight','normal')
%    txt = sprintf('D_0 = %2.0f', D0);
%    text(6.2e4,800,txt,'fontsize',12)
%    txt = sprintf('k = %2.1e',k0);
%    text(6.2e4,350,txt,'fontsize',12)
%    set(gca,'fontsize',12)     
%      
% 
% figure(3)  % ----------------------------------------------------------
%    set(gcf,'units','normalized');
%    set(gcf,'position',[0.670 0.05 0.2 0.50]);
%    set(gcf,'color','w');
% 
% subplot(2,1,1)   
%    yP = Rd+Dd; xP = Idtot;
%    plot( xP,yP,'b+')
%    hold on
%    yP = Rm; xP = Itot;
%    plot( xP,yP,'r','linewidth',2)
%    grid on
%    xlabel('infections   I_{tot}')
%    ylabel('removals   R_m')
%    set(gca,'fontsize',12)   
% 
%  subplot(2,1,2)   
%    yP = Rd+Dd; xP = Id;
%    plot( xP,yP,'b+')
%    hold on
%    yP = Rm; xP = I;
%    plot( xP,yP,'r','linewidth',2)
%    grid on
%    xlabel('active infections   I')
%    ylabel('removals   R_m')
%    set(gca,'fontsize',12)   
%    
figure(9)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.4 0.1 0.40 0.40]);
  set(gcf,'color','w');

  subplot('position',[0.1 0.4 0.3 0.3])
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
   text(80000,1200,txt,'fontsize',12)
   txt = sprintf('k = %2.2e',k0);
   text(80000,800,txt,'fontsize',12)
   set(gca,'fontsize',12)
   
 subplot('position',[0.55 0.15 0.35 0.3])  
   yP = Rd+Dd; xP = Id;
   plot( xP,yP,'b+')
   hold on
   yP = Rm; xP = I;
   plot( xP,yP,'r','linewidth',2)
   grid on
   xlabel('active infections   I')
   ylabel('removals   R_m')
   set(gca,'fontsize',12)   

 subplot('position',[0.55 0.65 0.35 0.3])  
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
 
T = [28 28 0];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [44 44 0];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [57 57 0];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [73 73 0];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [85 85 0];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [110 109 1];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [197 194 3];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [284 279 5];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [394 385 5];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [473 464 5];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [598 581 6];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [806 786 9];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [1023 1000 12];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1155 1129 15];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1588 1557 20];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [1950 1913 26];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2329 2288 30];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [2808 2701 38];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2906 2758 41];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [3666 3503 56];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [4068 3901 60];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [4823 4639 77];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [5658 5311 97];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [6359 5699 111];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [7045 6237 133];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8088 7179 151];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [8939 8014 167];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10065 9100 195];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [11426 9978 222];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [12186 10353 248];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13205 11321 267];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [13640 11348 278];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14277 11968 295];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [15013 12088 345];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [16009 12484 375];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [16638 12557 404];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [17760 13131 439];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [18678 13404 469];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [19317 14012 499];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [19822 14511 505];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [20596 15262 528];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21458 15508 550];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [22393 15617 576];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [23170 14544 601];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24195 13568 641 ];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [28455 14341 654];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25516 13674 672];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [26171 14311 690];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27566 15031 749];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [28455 15146 802];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [29893 16537 849];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [31142 15377 874];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [31993 15561 888];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [33027 15491 915];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [34238 16613 960];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [35330 16637 1006];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [36550 16814 1029];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [37727 17044 1079];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [38642 13942 1111];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [39890 15168 1133];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [40855 16110 1153];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [41866 16436 1179];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [43502 18602 1217];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [44775 19002 1258];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [46787 19269 1308];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [47672 19462 1340 ];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [48677 19747 1360];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [49684 20307 1369];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [50672 20807 1402];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [51651 19852 1425];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [53507 20236 1468];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [54616 20827 1512];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [55485 21491 1527];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [56166 22138 1535];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [56693 19802 1542];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [57744 19726 1563];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [59121 19579 1602];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [61182 20647 1627];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [62705 20898 1670];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [63526 20196 1679];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [65729 21613 1693];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [67281 21047 1717];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [68134 21881 1736];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [70901 22305 1797];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [71469 22810 1802];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [74368 24660 1843];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [74928 24184 1849];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [76762 25142 1862];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [78542 25510 1892];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [81434 27041 1920];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [82786 28384 1930];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [84363 29938 1953];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [87275 28745 1985];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [88048 29521 1992];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [90898 29798 2011];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [93047 31898 2044];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [97123 34397 2082];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [101788 35837 2139];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [103613 37651 2147];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [110107 40819 2192];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [112142 42853 2193];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [118300 46894 2216];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [123333 50376 2243];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [126363 53326 2263];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [134867 58049 2322];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [138560 61417 2347];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [146792 66141 2403];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [148845 68191 2406];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [158157 74389 2433];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [160410 76615 2460];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [167269 79955 2496];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [178755 85472 2563];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [183044 89739 2585];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [193766 93698 2638];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [195769 95691 2648];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [206546 100073 2691];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [212985 101732 2768];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [219420 108112 2823];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [235923 114597 3003];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [241013 119641 3046];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [257294 126212 3202];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [259465 128357 3228];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [271442 131712 3311];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [283718 137882 3438];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [285776 139907 3471];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;% T = [];
T = [305162 145511 3714];   
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;% T = [];
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
       m1 = [2020 04 20];    % day zero
       m2 = [2020 05  1];    % 
       m3 = [2020 06  1];    % 
       m4 = [2020 07  1];    % 
       m5 = [2020 08  1];    % 
       m6 = [2020 09  1];    % 
       m7 = [2020 10  1];    % 
       m8 = [2020 11  1];    % 
   
      tm(1) = datenum(m1);
      tm(2) = datenum(m2) - tm(1);
      tm(3) = datenum(m3) - tm(1);
      tm(4) = datenum(m4) - tm(1);
      tm(5) = datenum(m5) - tm(1);
      tm(6) = datenum(m6) - tm(1);
      tm(7) = datenum(m7) - tm(1);
      tm(8) = datenum(m8) - tm(1);
   
%       if z == 1; yMax = Imax; end
%       if z == 2; yMax = Rmax; end
%       if z == 3; yMax = Dmax; end
%       if z == 4; yMax = ITmax; end
%       if z == 5; yMax = Smax; end
      
      for n = 2:8
        plot([tm(n), tm(n)], [0, yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
        hold on
      end
      
      ys = 0.9;
      M = 'AMJJAS0';
      for c = 1:6
        text(5+tm(c+1), ys*yMax,M(c))
      end

end
