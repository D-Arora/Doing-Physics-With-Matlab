% covid19CHINA.m

% Revised:  2020 05 16

% Mathemtatical model for the spread of the covid19 virus
% Model for AUSTRALIAL spread of virus
% Data from https://www.worldometers.info/coronavirus/
% Date file:  covid19.mat  variable covid19
% Model specified by the variables in the INPUT SECTION of the Script

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/sird19EA.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b

% Populations: Data d / Model
% Susceptible individuals:  Sd  S
% Active infections (currently infected);  Id  I
% Removed Rm
% Recovered:  Rd  R
% Dead:  Dd  D


close all
clear
clc

global Imax Rmax Dmax ITmax Smax

% DATA  ===============================================================
load('covid19.mat')
Ndays  = find(covid19 == 0,1) - 1;  

% save covid19.mat covid19

% DATA: variable covid19
% Col 1: World   Id   Currently Infected Cases
% Col 2: World   Rd   Recovered
% Col 3: World   Dd   Deaths

% Data for China with Script svr19CHINA.m
% Col 4: China   Total infected cases
% Col 5: China   Currently Infected Cases
% Col 6: China   Deaths
% Col 7: USA     Total infected cases
% Col 8: USA     Currently infected cases
% Col 9: USA     Deaths
% Col 10: AUST     Total infected cases
% Col 11: AUST     Currently infected cases
% Col 12: AUST     Deaths
% Col 13: IND     Total infected cases
% Col 14: IND     Currently infected cases
% Col 15: IND     Deaths

% DATA: Id Rd Dd / Number of data days for COVID19  
  Idtot  = covid19(53:Ndays,13);
  Id     = covid19(53:Ndays,14);
  Dd     = covid19(53:Ndays,15);  
  
  Ndays = length(Idtot);
  
  Rd = Idtot - Id - Dd;
  Idmax = max(Id);
  Rdmax = max(Rd);
  Ddmax = max(Dd);
  ITdmax = max(Idtot);

% MODEL INPUT PARAMETERS  ============================================
% Model: number of time steps 
  num = 5000;
% Maximum number of days for model   [200]
  tMax = 200;

% Initial number of infections   [0.004] 
  I = zeros(num,1);
  I(1) = 1e-4;
% Adjustable model parameters  
% S --> I  [0.38] 
  a = 0.20;
% I --> Rm  [0.036]
  b = 0.045;
%  r = 0.94;    % recoveries
%  d = 1 - r;   % deaths
% Population scaling factor
%  f = 8.7171e4;
  f = 1.5e5;
% Time index
  t_index = 77; 
% Exponential growth
     th = 5; t0 = -6; 
% Starting Date
  zSTART = datenum([2020 3 14]);

% SETUP  ==============================================================
% Time [days]
  t = linspace(0,tMax,num)';
  dt = t(2) - t(1);
% Exponential growth
  k = log(2)/th;  tE = 0:80;
  v0 = 177400/exp(k*60);
  v = v0.*exp(k*(tE - t0));

% Model: S  D  R  
  S  = zeros(num,1);
  Rm = zeros(num,1);
    
% Model: Initial values for populations
  S(1) = 1;
  
% MODEL: Time Evolution  ==============================================
for n = 1:num-1
    
   if n == 400; S(n)  = 0.90; end
    if n == 1600; S(n) = 0.60; end
    if n == 1750; S(n) = 0.65; end
    if n == 1900; S(n) = 0.50; end
    if n == 2000; S(n) = 0.50; end
    if n == 2100; S(n) = 0.6; end 
    if n == 2200; S(n) = 0.6; end 
    if n == 2300; S(n) = 0.6; end 
    
  dS = -a*S(n)*I(n)*dt;
  dI = -dS;
  dRm = b*I(n)*dt;
  S(n+1) =  S(n) + dS;
   
%   if n == 400; S(n+1) = 1.0; end
%   if n == 1600; S(n+1) = 0.70; end
%   if n == 1750; S(n+1) = 0.60; end
%   if n == 1900; S(n+1) = 0.50; end
%   if n == 2000; S(n+1) = 0.50; end
%   if n == 2100; S(n+1) = 0.50; end
  
  I(n+1) =  I(n) + dI - dRm;
  Rm(n+1) = Rm(n)+ dRm;
end
 

  I = f.*I;
  Rm = f.*Rm;
  k0 = 1.4e-5;
  D0 = 11e3;
  D = D0.*(1 - exp(-k0.*Rm));
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
  
% Percentages: Data
  pId = 100*Id(end) / Idtot(end);
  pRd = 100*Rd(end) / Idtot(end);
  pDd = 100*Dd(end) / Idtot(end);
  

% *********************************************************************  
% GRAPHICS
% ********************************************************************* 

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.1 0.71 0.6]);
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
    
    xlim([0 tMax])
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
  %  set(gca,'Position', [0.5703-0.04 0.7093+0.04 0.3347 0.2157]);
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
   
    xlim([0 tMax])
    grid on
    ylabel('R')
    xlabel('days elapsed')
    title('Recoveries  R','fontweight','normal')
    set(gca,'fontsize',12)
 %   set(gca,'Position', [0.0800 0.4096+0.04 0.3347 0.2157]);
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
    
    xlim([0 tMax])
    grid on
    ylabel('D')
    xlabel('days elapsed')
    title('Deaths  D','fontweight','normal')
    set(gca,'fontsize',12) 
 %   set(gca,'Position', [0.5703-0.04 0.4096+0.04 0.3347 0.2157]);
    yMax = ylim;
    z = yMax(2);
    plotMonths(z)
 

 subplot(2,3,5)
    xP = t; yP = S;
    plot(xP,yP,'b','linewidth',2)
    hold on
   
    xlim([0 tMax])
    ylim([0 1.1*Smax])
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
   
   txt = 'Start date    14 March 2020';
   Htext = text(0,hh,txt,'fontsize',12);
   set(Htext,'color','k')
   
   hh = hh+2*dh;
   txt = 'MODEL PARAMETERS';
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('I(1) = %3.0e   f = %2.2e   a = %2.3f   b = %2.3f', I(1)/f, f,a ,b);
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
   
     hh = hh+dh; 
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
%     set(gca,'Position', [0.5703-0.08 0.1100+0.04 0.3347 0.2157]);
     
     
%  figure(2)  % ----------------------------------------------------------
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
%    text(1.52e5,3e3,txt,'fontsize',12)
%    txt = sprintf('k = %2.1e',k0);
%    text(1.52e5,1e3,txt,'fontsize',12)
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
   text(2.2e5,4000,txt,'fontsize',12)
   txt = sprintf('k = %2.2e',k0);
   text(2.2e5,2000,txt,'fontsize',12)
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

      
   
   
% =====================================================================
     
function  plotMonths(z)
     
       yMax = z;
     % Months
       m1 = [2020 03 14];    % day zero
       m2 = [2020 04  1];    % Apr
       m3 = [2020 05  1];    % May
       m4 = [2020 06  1];    % Jun
       m5 = [2020 07  1];    % Jul
       m6 = [2020 08  1];    % Aug
       m7 = [2020 09  1];    % Sep
       m8 = [2020 10  1];    % Oct
   
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
      text(22, ys*yMax,'A')
      text(52, ys*yMax,'M')
      text(82, ys*yMax,'J')
      text(112,ys*yMax,'J')
      text(142,ys*yMax,'A')
      text(172,ys*yMax,'S')
      
   end
   