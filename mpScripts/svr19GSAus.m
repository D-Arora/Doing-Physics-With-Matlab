% svr19.m

% Revised:  2020 03 25

% Mathemtatical model for the spread of the covid19 virus
% Model for GLOBAL spread of virus
% Model for CHINA spread of virus: download svr19CHINA.m
% Data from https://www.worldometers.info/coronavirus/
% Date file:  covid19.mat  variable covid19
% Model specified by the variables in the INPUT SECTION of the Script

% DOING PHYSICS WITH MATLAB: 
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/sird19E.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% Matlab 2019b

% Populations: Data d / Model
% Susceptible individuals:  Sd  S
% Active infections (currently infected);  Id  I
% Recovered:  Rd  R
% Dead:  Dd  D

close all
clear
clc

global Imax Rmax Dmax ITmax Smax

load('covid19.mat')

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


% Initialize arrays  **************************************************
  % Model: number of time steps 
  num = 5000;
  I = zeros(num,1);
  N = zeros(num,1);
  D = zeros(num,1);
  R = zeros(num,1);

% INPUTS: Model Parameters ********************************************
% Maximum number of days for model   [200]
  tMax = 100;
% Adjustible model parameters  
  a(1) = 0.35;
  a(2) = 0.042;
  a(3) = 0.985;
  a(4) = 1;
  I(1) = 100;
  S0 = 1e5;
  w  = 11.55;
  t0 = 0;
  
  
% SETUP  ************************************************************** 
% Time [days]
  t = linspace(0,tMax,num)';
  dt = t(2) - t(1);

% Susceptible Population
  S = S0 .* exp(-0.5.*((t-t0)./w).^2);

% Exponential increase
  th = 3.5; k = log(2)/th;
  Iexp = I(1).*exp(k*(t));  
  
  
% DATA: Id Rd Dd / Number of data days for COVID19 ******************** 
  Ndays  = find(covid19 == 0,1) - 1;  
  DayS = 50;
  Idtot  = covid19(DayS:Ndays,10);
  Id     = covid19(DayS:Ndays,11);
  Dd     = covid19(DayS:Ndays,12);
  Ld     = length(Id);

  Rd = Idtot - Id - Dd;
  Idmax = max(Id);
  Rdmax = max(Rd);
  Ddmax = max(Dd);
  ITdmax = max(Idtot);
  
% DATA: Percentages
  pId = 100*Id(end) / Idtot(end);
  pRd = 100*Rd(end) / Idtot(end);
  pDd = 100*Dd(end) / Idtot(end);


% MODEL  **************************************************************  
for n = 1:num-1
  dS = a(1)*(S(n)/S0)*I(n)*dt;
  dN = a(2)*I(n)*dt;
  I(n+1) = I(n) + dS - dN;
  N(n+1) = N(n) + dN;
  R(n+1) = a(3)*N(n+1);
  D(n+1) = N(n+1) - R(n+1);
 % S(n+1) = S(n+1) - dS;%a(4)*dS;
end

  Itot = I + R + D; 
  Imax = max(I);
  ITmax= max(Itot);
  Rmax = max(R);
  Dmax = max(D);
  Smax = max(S);
  
% Percentages: model
  pI = 100*I(end) / Itot(end);
  pR = 100*R(end) / Itot(end);
  pD = 100*D(end) / Itot(end);  
  
 
% GRAPHICS  *********************************************************** 

figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.3 0.75]);
  set(gcf,'color','w');

subplot(4,1,1)  % INFECTIONS
    xP = 1: Ld;
    yP = Id;
    plot(xP,yP,'b+')   

    hold on
    
    xP = t; yP = I;
    plot(xP,yP,'b','linewidth',2)
   
     plotMonths(1)
     
    grid on
    title('Active Infections  I')
    ylabel('I','FontName','Times New Roman')  
    set(gca,'fontsize',12)
    xlim([0 tMax])
    
subplot(4,1,2)   % RECOVERIES
    xP = 1: Ld;
    yP = Rd; 
    plot(xP,yP,'k+')
    hold on
 
    xP = t; yP = R;
    plot(xP,yP,'k','linewidth',2)
    
    plotMonths(2)
    
    grid on

    ylabel('R')
    title('Recoveries  R')
    set(gca,'fontsize',12)
    xlim([0 tMax])
 
 subplot(4,1,3)    % DEATHS
    xP = 1: Ld;
    yP = Dd; 
    plot(xP,yP,'r+')
    hold on
 
    xP = t; yP = D;
    plot(xP,yP,'r','linewidth',2)
    
    plotMonths(3)

    grid on
    ylabel('D')
    title('Deaths  D')
    set(gca,'fontsize',12) 
    xlim([0 tMax])
    
 subplot(4,1,4)   % Total infections
    xP = t(1:2500);
    yP = Iexp(1:2500);
    plot(xP,yP,'r','linewidth',1.2)
     ylim([0 1.2*ITmax])
    hold on
 
    xP = 1: Ld;
    yP = Idtot; 
    plot(xP,yP,'b+')
    hold on
    xP = t; yP = Itot;
    plot(xP,yP,'b','linewidth',2)
    
    plotMonths(4)
    
    grid on
    xlabel('days elapsed')
    ylabel('I_{tot}','FontName','Times New Roman')
    title('Total Infections I_{tot}');
    set(gca,'fontsize',12) 
    xlim([0 tMax])

figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.60 0.05 0.25 0.25]);
   set(gcf,'color','w');

   plotMonths(5)   

   xP = t; yP = S;
   plot(xP,yP,'k','linewidth',2)
   grid on
   xlabel('days elapsed')
   ylabel('Susceptible population  S') 
      
   txt = sprintf('S_{max} = %3.0e',Smax);
   text(100,11.2e4,txt,'fontsize',12)

    
    
    
figure(7)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.33 0.1 0.25 0.5]);
   set(gcf,'color','w');

   xlim([0 120])
   ylim([-60 250])
   hh = 250; dh = -25;
   
   txt = 'Start date:  22 January 2020';
   Htext = text(0,hh,txt,'fontsize',12);
   set(Htext,'color','k')
   
   hh = hh+dh;
   txt = 'Model Parameters';
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('I(1) = %3.0f   S_0 = %2.0e   w = %2.2f   t_0 = %2.0f', ...
       I(1),S0,w,t0);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('a(1) = %3.3f   a(2) = %2.3f   a(3) = %2.3f   a(4)= %2.3f', ...
       a(1),a(2),a(3),a(4));
   text(0,hh,txt,'fontsize',12)
      
   hh = hh+2*dh; 
   txt = 'DATA';
   text(0,hh,txt,'fontsize',12)
   
    z = datenum([2020 1 22]);
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
    pR = 100*Rd(end)/(Rd(end)+Dd(end));
    pD = 100*Dd(end)/(Rd(end)+Dd(end));
    hh = hh+dh;
    txt = sprintf('percent: Active = %2.1f  Recoveries = %2.1f    Deaths = %2.1f   ',pId,pRd,pDd);
    text(2,hh,txt,'fontsize',12)  
   
   hh = hh+2*dh; 
   txt = 'MODEL PREDICTIONS';
   text(0,hh,txt,'fontsize',12)
   
    z = datenum([2020 1 22]);
    z = z + tMax;
    z = datetime(z,'ConvertFrom','datenum');
    z.Format = 'dd-MMM-yyyy';
    txt = cellstr(z) ; 
    text(60,hh,txt,'fontsize',12)
   
   hh = hh+dh;
    z1 = I(end); z2 = R(end); z3 = D(end); z4 = z1+z2+z3;  
    txt = sprintf('I = %5.0f    R = %5.0f    D = %5.0f    I_{tot} = %5.0f', z1,z2,z3,z4);
    text(6,hh,txt,'fontsize',12)

   hh = hh+dh;
    txt = 'Peak';
    text(10,hh,txt,'fontsize',12)
    z = find(I == Imax,1);
    z = t(z) + datenum([2020 1 22]);
    z = datetime(z,'ConvertFrom','datenum');
    z.Format = 'dd-MMM-yyyy';
    txt = cellstr(z) ; 
    text(30,hh,txt,'fontsize',12)
   
    txt = sprintf('I_{peak} = %5.0f',Imax);
    text(70,hh,txt,'fontsize',12)
   
    hh = hh+dh;
    txt = sprintf('percent: Active = %2.1f  Recoveries = %2.1f    Deaths = %2.1f   ',pI,pR,pD);
    text(2,hh,txt,'fontsize',12)   
       

    
    axis off
   
    
    function  pM = plotMonths(z)
     global Imax Rmax Dmax ITmax Smax
     
     % Months
       m1 = [2020 01 22];    % day zero
       m2 = [2020 02  1];    % Feb
       m3 = [2020 03  1];    % Mar
       m4 = [2020 04  1];    % Apr
       m5 = [2020 05  1];    % May
       m6 = [2020 06  1];    % Jun
       m7 = [2020 07  1];    % Jul
       m8 = [2020 08  1];    % Aug
   
      tm(1) = datenum(m1);
      tm(2) = datenum(m2) - tm(1);
      tm(3) = datenum(m3) - tm(1);
      tm(4) = datenum(m4) - tm(1);
      tm(5) = datenum(m5) - tm(1);
      tm(6) = datenum(m6) - tm(1);
      tm(7) = datenum(m7) - tm(1);
      tm(8) = datenum(m8) - tm(1);
   
      if z == 1; yMax = Imax; end
      if z == 2; yMax = Rmax; end
      if z == 3; yMax = Dmax; end
      if z == 4; yMax = ITmax; end
      if z == 5; yMax = Smax; end
      
      for n = 2:8
        plot([tm(n), tm(n)], [0, yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
        hold on
      end
      
      ys = 0.8;
      text(22, ys*yMax,'F')
      text(52, ys*yMax,'M')
      text(82, ys*yMax,'A')
      text(112,ys*yMax,'M')
      text(142,ys*yMax,'J')
      text(172,ys*yMax,'J')
      
   end
   
   