% svr19.m

% Revised:  2020 03 25

% Mathemtatical model for the spread of the covid19 virus
% Model for GLOBAL spread of virus
% Model for CHINA spread of virus: download svr19CHINA.m
% Data from https://www.worldometers.info/coronavirus/
% Date file:  covid19.mat  variable covid19
% Model specified by the variables in the INPUT SECTION of the Script

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/sird19E.htm
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


close all
clear
clc

global Imax Rmax Dmax ITmax Smax f

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


% Model: number of time steps 
  num = 5000;
% Maximum number of days for model   [200]
  tMax = 200;
% TIME [days]
  t = linspace(0,tMax,num)';
  dt = t(2) - t(1);
    
  flagS = input('Select: 1 World / 2 China / 3 USA / 4 Australia  ');

  switch flagS    % WORLD
  case 1
     %  Initial number of infections   [0.004] 
        I = zeros(num,1);
        I(1) = 0.004;
     % Adjustible model parameters  
     % S --> I  [0.38] 
        a = 0.39;
     % I --> Rm  [0.036]
        b = 0.028;
        r = 0.85;     % recoveries
        d = 1 - r;    % deaths
     % Time index
       t_index = 90;%27;
     % Exponential growth
       th = 5; t0 = 0;
     
     % SURGE FACTORS
       s = zeros(num,1);
       for n = 1:num-1
          if t(n) > 30;  s(n) = 0.33; end
          if t(n) > 58;  s(n) = 0.85; end 
          if t(n) > 64;  s(n) = 1.90; end   
          if t(n) > 68;  s(n) = 2.80; end 
          if t(n) > 71;  s(n) = 3.50; end  
          if t(n) > 75;  s(n) = 4.50; end
          if t(n) > 78;  s(n) = 5.50; end
          if t(n) > 80;  s(n) = 6.50; end
          if t(n) > 85;  s(n) = 7.50; end
          if t(n) > 90;  s(n) = 8.50; end
          if t(n) > 95;  s(n) = 9.00; end
          if t(n) > 98; s(n) = 10.00; end
          if t(n) > 106; s(n) = 10.5; end
          if t(n) > 109; s(n) = 11; end
          if t(n) > 120; s(n) = 0.00; end  
       end    
             
     % DATA: Id Rd Dd / Number of data days for COVID19  
       Id = covid19(1:Ndays,1);
       Rd = covid19(1:Ndays,2);
       Dd = covid19(1:Ndays,3);
     
       Idtot = Id + Rd + Dd;
       Idmax = max(Id);
       Rdmax = max(Rd);
       Ddmax = max(Dd);
       ITdmax = max(Idtot);
    
       
  case 2    % CHINA
      % Initial number of infections   [0.004] 
        I = zeros(num,1);
        I(1) = 0.004;
      % Adjustible model parameters  
      
        a = 0.35;    % S --> I  [0.35] 
      % I --> Rm  [0.036]
        b = 0.035;
        r = 0.94;     % recoveries
        d = 1 - r;    % deaths
      % Time index
         t_index = 27;
      % Exponential growth
       th = 5; t0 = -30;
       
       % SURGE FACTORS
        s = zeros(num,1);
        for n = 1:num-1
           if t(n) > 90 && t(n) < 120
           s(n) = 0.0;
           end
        end  
        
      % DATA: Id Rd Dd / Number of data days for COVID19 
         Idtot  = covid19(1:Ndays,4);
         Id     = covid19(1:Ndays,5);
         Dd     = covid19(1:Ndays,6);
     
         Rd = Idtot - Id - Dd;
         Idmax = max(Id);
         Rdmax = max(Rd);
         Ddmax = max(Dd);
         ITdmax = max(Idtot);
 
 case 3     % USA      
      % Initial number of infections   [0.004] 
        I = zeros(num,1);
        I(1) = 0.001;
      % Adjustible model parameters  
      % S --> I  [0.38] 
        a = 0.1;
      % I --> Rm  [0.036]
        b = 0.012;
        r = 0.80;     % recoveries
        d = 1 - r;    % deaths
      % Time index
        t_index = 57;
      % Exponential growth
        th = 5; t0 = 10;   
        
   
      % SURGE FACTORS
       s = zeros(num,1);
         for n = 1:num-1
            if t(n) > 50; s(n) = 0.19; end
            if t(n) > 62; s(n) = 0.35; end 
            if t(n) > 70; s(n) = 1.50; end 
            if t(n) > 77; s(n) = 2.20; end
            if t(n) > 81; s(n) = 3.30; end
            if t(n) > 87; s(n) = 4.50; end
            if t(n) > 94; s(n) = 5.20; end  
            if t(n) > 100; s(n) = 6.00; end 
            if t(n) > 105; s(n) = 6.20; end 
            if t(n) > 120; s(n) = 0.00; end  
         end 
         
         
      % DATA: Id Rd Dd / Number of data days for COVID19  
        Idtot  = covid19(1:Ndays,7);
        Id     = covid19(1:Ndays,8);
        Dd     = covid19(1:Ndays,9);  
        
        Rd = Idtot - Id - Dd;
        Idmax = max(Id);
        Rdmax = max(Rd);
        Ddmax = max(Dd);
        ITdmax = max(Idtot);  
     
         
  case 4    % AUSTRALIA
    % Initial number of infections   [0.004] 
      I = zeros(num,1);
      I(1) = 0.004;
   % Adjustable model parameters  
   % S --> I  [0.38] 
     a = 0.1;
   % I --> Rm  [0.036]
     b = 0.07;
     r = 0.984;    % recoveries
     d = 1 - r;   % deaths
   % Time index
     t_index = 77; 
   % Exponential growth
     th = 5; t0 = 35; 
   
  % SURGE FACTORS
    s = zeros(num,1);
    for n = 1:num-1
      if t(n) > 47; s(n) = 0.14; end
      if t(n) > 58; s(n) = 0   ; end
    end   
    
   % DATA: Id Rd Dd / Number of data days for COVID19  
     Idtot  = covid19(1:Ndays,10);
     Id     = covid19(1:Ndays,11);
     Dd     = covid19(1:Ndays,12);  
     
     Rd = Idtot - Id - Dd;
     Idmax = max(Id);
     Rdmax = max(Rd);
     Ddmax = max(Dd);
     ITdmax = max(Idtot);

  end


% *********************************************************************  
% SETUP
% ********************************************************************* 

% Exponential growth
  k = log(2)/th;
  v0 = 177400/exp(k*60);
  v = v0.*exp(k*(t - t0));

% Model: S  D  R  
  S  = zeros(num,1);
  Rm = zeros(num,1);
    
% Model: Initial values for populations
  S(1) = 1;
  
% *********************************************************************  
% MODEL: Time Evolution
% ********************************************************************* 
 
for n = 1:num-1
  dS = -a*S(n)*I(n)*dt;
  dI = -dS;
  dRm = b*I(n)*dt;
  surge = s(n)*S(n)*dt;
  S(n+1) =  S(n) + dS + surge;
  I(n+1) =  I(n) + dI - dRm;
  Rm(n+1) = Rm(n)+ dRm;
end
  
  R = r.*Rm;    % recoveries
  D = d.*Rm;    % deaths
  Itot = I + R + D;
  
% Scaling I, R and D values at day specified by t_index 
    index = find(t > t_index,1);
    f = Id(t_index)./ I(index);
%   f = 1.6e3;
% Scale populations
  I = f.*I;
  R = f.*R;
  D = f*D;
  Itot = f.*Itot;
% 
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
  set(gcf,'position',[0.02 0.05 0.3 0.75]);
  set(gcf,'color','w');

subplot(4,1,1)  % INFECTIONS
    xP = 1: Ndays;
    yP = Id;
    plot(xP,yP,'b+')   
    hold on
    
    xP = t_index; yP = Id(t_index);
    plot(xP,yP,'bo','markerfacecolor','b','markersize',8)
    
    xP = t; yP = I;
    plot(xP,yP,'b','linewidth',2)
   
    plotMonths(1)
    
    grid on
    title('Active Infections  I')
    ylabel('I','FontName','Times New Roman')  
    set(gca,'fontsize',12)
  
    
subplot(4,1,2)   % RECOVERIES
    xP = 1: Ndays;
    yP = Rd; 
    plot(xP,yP,'k+')
    hold on
 
    xP = t_index; yP = Rd(t_index);
    plot(xP,yP,'ko','markerfacecolor','k','markersize',8)
    
    xP = t; yP = R;
    plot(xP,yP,'k','linewidth',2)
    
    plotMonths(2)
    
    grid on

    ylabel('R')
    title('Recoveries  R')
    set(gca,'fontsize',12)
   
 
 subplot(4,1,3)    % DEATHS
    xP = 1: Ndays;
    yP = Dd; 
    plot(xP,yP,'r+')
    hold on
 
    xP = t_index; yP = Dd(t_index);
    plot(xP,yP,'ro','markerfacecolor','r','markersize',8)
    
    xP = t; yP = D;
    plot(xP,yP,'r','linewidth',2)
    
    plotMonths(3)
    grid on
    
    ylabel('D')
    title('Deaths  D')
    set(gca,'fontsize',12) 
    
    
 subplot(4,1,4)   % Total infections
    xP = t(1:2500);
    yP = v(1:2500);
    plot(xP,yP,'r','linewidth',1.2)
   
    hold on
    xP = 1: Ndays;
    yP = Idtot; 
    plot(xP,yP,'b+')
    hold on
    xP = t; yP = Itot;
    plot(xP,yP,'b','linewidth',2)
    
    plotMonths(4)
    
    ylim([0 1.2*ITmax])
    grid on
    xlabel('days elapsed')
    ylabel('I_{tot}','FontName','Times New Roman')
    title('Total Infections I_{tot}');
    set(gca,'fontsize',12) 
        
    
figure(2)  % ----------------------------------------------------------
   set(gcf,'units','normalized');
   set(gcf,'position',[0.60 0.05 0.25 0.25]);
   set(gcf,'color','w');

   fig = figure(2);
   left_color = [0 0 0];
   right_color = [1 0 1];
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
   
    plotMonths(5)   
        
    yyaxis left
    xP = t; yP = f.*S;
    plot(xP,yP,'k','linewidth',2)
    grid on
    xlabel('days elapsed')
    ylabel('Susceptible population  S') 
    
    %ylim([0 1])
    yyaxis right
    xP = t; yP = s;
    plot(xP,yP,'m','linewidth',1)
    grid on
    ylabel('surge factor  s') 
        
    set(gca,'fontsize',12)    
    
    
figure(7)  % ----------------------------------------------------------
   set(gcf,'units','normalized');
   set(gcf,'position',[0.33 0.05 0.27 0.35]);
   set(gcf,'color','w');

   xlim([0 120])
   ylim([0 250])
   hh = 250; dh = -25;
   
   txt = 'Start date:  22 January 2020';
   Htext = text(0,hh,txt,'fontsize',12);
   set(Htext,'color','k')
   
   hh = hh+dh;
   txt = 'Model Parameters';
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = sprintf('I(0) = %3.0f   a = %2.3f   b = %2.3f   d = %2.3f', ...
       I(1),a,b,d);
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
   
% =====================================================================
     
function  pM = plotMonths(z)
     global Imax Rmax Dmax ITmax Smax f
     
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
      if z == 5; yMax = f.*Smax; end
      
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
   
   