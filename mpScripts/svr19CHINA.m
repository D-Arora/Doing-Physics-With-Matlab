% svr19CHINA.m

% Revised:  2020 03 25

% Mathemtatical model for the spread of the covid19 virus
% Model for GLOBAL spread of virus
% Model for CHINA spread of virus: download svr19CHINA.m
% Data from https://www.worldometers.info/coronavirus/
% Date file:  covid19.mat  variable covid19
% Model specified by the variables in the INPUT SECTION of the Script

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
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

% Model: number of time steps 
  num = 5000;

% *********************************************************************  
% INPUTS: Model Parameters
% ********************************************************************* 

% Maximum number of days for model   [200]
  tMax = 200;
% Initial number of infections   [0.004] 
  I = zeros(num,1);
  I(1) = 0.004;
% Adjustible model parameters  
% S --> I  [0.38] 
  a = 0.35;
% I --> R  [0.036]
  b = 0.038;
% I --> D  [0.0025]
  d = 0.0015;

% See SETUP for input of surge factors

% *********************************************************************  
% SETUP
% ********************************************************************* 
  
% DATA: Id Rd Dd / Number of data days for COVID19  
  Ndays  = find(covid19 == 0,1) - 1;  
  Idtot  = covid19(1:Ndays,4);
  Id     = covid19(1:Ndays,5);
  Dd     = covid19(1:Ndays,6);

% Model: S  D  R  
  S = zeros(num,1);
  D = zeros(num,1);
  R = zeros(num,1);
  
% Model: Initial values for populations
  S(1) = 1;
  D(1) = 0;
  R(1) = 0;   
  
% TIME [days]
  t = linspace(0,tMax,num)';
  dt = t(2) - t(1);

  
% SURGE FACTORS
  s = zeros(num,1);
   for n = 1:num-1
     if t(n) > 90 && t(n) < 120
        s(n) = 0.0;
     end
   end

   
% *********************************************************************  
% MODEL: Time Evolution
% ********************************************************************* 
 
for n = 1:num-1
  dS = -a*S(n)*I(n)*dt;
  dI = -dS;
  I(n+1) = I(n) + dI;
  dR = b*I(n)*dt;
  dD = d*I(n)*dt;
  surge = s(n)*S(n)*dt;
  R(n+1) = R(n)+ dR;
  D(n+1) = D(n) + dD;
  S(n+1) = S(n) + dS + surge;
  I(n+1) = I(n+1) - dR - dD;
end

% Scaling I, R and D values at day specified by t_index 
%   I(27) = Id(27)
  t_index = 27;
  index = find(t > t_index,1);
  f = Id(t_index)./ I(index);

% Scale populations
  I = f.*I;
  R = f.*R;
  D = f*D;
  
  Imax = max(I);
  Rmax = max(R);
  Dmax = max(D);
  Smax = max(S);
  
% Model: total population that has been infected;  
  Itot = I + R + D;
  ITmax = max(Itot);
  
% Data
  Rd = Idtot - Id - Dd;
  Idmax = max(Id);
  Rdmax = max(Rd);
  Ddmax = max(Dd);
  ITdmax = max(Idtot);
  
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
  
    hold on    
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
    xP = 1: Ndays;
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
    

figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.60 0.05 0.25 0.25]);
   set(gcf,'color','w');

   fig = figure(2);
   left_color = [0 0 0];
   right_color = [1 0 1];
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
   
    hold on
    plotMonths(5)   
        
    yyaxis left
    xP = t; yP = S;
    plot(xP,yP,'k','linewidth',2)
    grid on
    xlabel('days elapsed')
    ylabel('Susceptible population  S') 
      
    yyaxis right
    xP = t; yP = s;
    plot(xP,yP,'m','linewidth',1)
    grid on
    ylabel('surge factor  s') 
        
    set(gca,'fontsize',12)    
    box on
    
figure(7)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.33 0.05 0.25 0.35]);
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
   
      if z == 1; yMax = Imax;  end
      if z == 2; yMax = Rmax;  end
      if z == 3; yMax = Dmax;  end
      if z == 4; yMax = ITmax; end
      if z == 5; yMax = Smax;  end
      
      for n = 2:8
        plot([tm(n), tm(n)], [0, yMax],'color', [0.5 0.5 0.5],'linewidth',1);
      end
      
      ys = 0.8;
      text(22, ys*yMax,'F')
      text(52, ys*yMax,'M')
      text(82, ys*yMax,'A')
      text(112,ys*yMax,'M')
      text(142,ys*yMax,'J')
      text(172,ys*yMax,'J')
      
   end