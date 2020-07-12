% covid19UK.m

clear 
clc
close all

global Imax Rmax Dmax ITmax Smax

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
  

% MODEL  ============================================================
% Model: number of time steps 
  num = 5000;
% Maximum number of days for model   [200]
  tMax = 200;
% TIME [days]
  t = linspace(0,tMax,num)';
  dt = t(2) - t(1);

  S  = zeros(num,1);
  Rm = zeros(num,1);
   
 % Initial number of infections   [0.004] 
   I = zeros(num,1);
   I(1) = 0.004;
   S(1) = 1;
 % Adjustible model parameters  
 % S --> I  [0.38] 
   a = 0.13;
 % I --> Rm  [0.036]
   b = 0.03;
   r = 0.85;     % recoveries
   d = 1 - r;    % deaths
   f = 4.5e5;
 % Time index
    t_index = 50;  
   
% SURGE FACTORS
       s = zeros(num,1);
%        for n = 1:num-1
%           if t(n) > 30;  s(n) = 0.33; end
%           if t(n) > 58;  s(n) = 0.85; end 
%           if t(n) > 64;  s(n) = 1.90; end   
%           if t(n) > 68;  s(n) = 2.80; end 
%           if t(n) > 71;  s(n) = 3.50; end  
%           if t(n) > 75;  s(n) = 4.50; end
%           if t(n) > 78;  s(n) = 6.50; end
%           if t(n) > 81;  s(n) = 7.50; end
%           if t(n) > 90; s(n) = 0.00; end  
%        end    
  
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
% index = find(t > t_index,1);
%  f = Id(t_index)./ I(index);
 
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

% Exponential growth
  th = 4.4; t0 = 0; tE = 0 : 40;%Ndays;
  k = log(2)/th;
  v0 = Idtot(1);
  v = v0.*exp(k*(tE - t0));    
  
% Percentages: model
  pI = 100*I(end) / Itot(end);
  pR = 100*R(end) / Itot(end);
  pD = 100*D(end) / Itot(end);
  
% Percentages: Data
  pId = 100*Id(end) / Idtot(end);
  pRd = 100*Rd(end) / Idtot(end);
  pDd = 100*Dd(end) / Idtot(end);
  
  
% GRAPHICS  ===========================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.05 0.3 0.75]);
  set(gcf,'color','w');
  
  subplot(4,1,1)  % INFECTIONS
    xP = 1: Ndays;
    yP = Id;
    plot(xP,yP,'b+')   
    hold on
    
    xP = t; yP = I;
    plot(xP,yP,'b','linewidth',2)
%    
    plotMonths(1)
%     
    grid on
    xlim([0 tMax])
    title('Active Infections  I')
    ylabel('I','FontName','Times New Roman')  
    set(gca,'fontsize',12)

subplot(4,1,2)   % RECOVERIES
    xP = 1: Ndays;
    yP = Rd; 
    plot(xP,yP,'k+')
    hold on

    xP = t; yP = R;
    plot(xP,yP,'k','linewidth',2)
    
    plotMonths(2)
    
    grid on
    xlim([0 tMax])
    ylabel('R')
    title('Recoveries  R')
    set(gca,'fontsize',12)
   
 
 subplot(4,1,3)    % DEATHS
    xP = 1: Ndays;
    yP = Dd; 
    plot(xP,yP,'r+')
    hold on
    
    xP = t; yP = D;
    plot(xP,yP,'r','linewidth',2)
     
     plotMonths(3)
    grid on
    xlim([0 tMax])
    ylabel('D')
    title('Deaths  D')
    set(gca,'fontsize',12) 
    
    
 subplot(4,1,4)   % Total infections
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
     
    plotMonths(4)

    grid on
    xlim([0 tMax])
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
      
    yyaxis right
    xP = t; yP = s;
    plot(xP,yP,'m','linewidth',1)
    grid on
    ylabel('surge factor  s') 
        
    set(gca,'fontsize',12)    
    xlim([0 tMax])
 
figure(7)  % ----------------------------------------------------------
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
   
   z = datenum([2020 3 12]);
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
   
    z = datenum([2020 3 12]);
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
     z = t(z) + datenum([2020 3 12]);
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
       
  
% FUNCTIONS ============================================================
  
function covid = sir
c = 1;
T = [590 562 10];
 covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [798 769 11];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1140 1101 21];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1391 1336 35];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1543 1436 55];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1950 1814 71];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2626 2457 104];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [3269 3060 144];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [3983 3741 177];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [5018 4692 233];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [5683 5309 281];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [6650 6180 335];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8077 7520 422];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9529 8931 463];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11658 10945 578];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14543 13649 759];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;     
T = [17089 15935 1019];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [19522 18159 1228];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;     
T = [22141 20598 1408];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;     
T = [25150 23226 1789];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [29474 26987 2352];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;     
T = [33718 30662 2921];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;  
T = [38168 34428 3605 ];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;      
T = [41903 37455 4313];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [47806 42737 4934];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [51608 46100 5373];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [55242 48948 6159];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;     
T = [60733 53501 7097];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [65077 56964 7978];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;     
T = [73758 64456 8958];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;     
T = [78991 68772 9875];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;     
T = [84279 73323 10612];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [88621 76948 11329];
     covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [93873 81422 12107];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [98476 85264 12868];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [108692 93772 14576];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [114217 98409 15464];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [120067 103663 16060];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [124743 107890 16509];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [129044 111363 17337];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [133495 115051 18100];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [138078 118996 18738];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [143464 123614 19506];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [148377 127714 20319];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [152840 131764 20732];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [157149 135713 21092];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [161145 139123 21678];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [165221 138780 26097];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [171253 144138 26771];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [177454 149600 27510];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [182260 153785 28131];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [186599 157809 28884];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [190584 161506 28734];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [194990 165129 29427];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [201101 170681 30076];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [206715 175756 30615];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [211364 179779 31241];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [215260 183329 31587];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [219183 186984 31855];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [223060 190651 32065];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [226463 193427 32692];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
%T = [];
%   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
%T = [];
%  
%T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
%T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
%T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
%T = [];
%  


end





function  plotMonths(z)
     global Imax Rmax Dmax ITmax Smax
     
     % Months
       m1 = [2020 03 12];    % day zero
       m2 = [2020 04  1];    % 
       m3 = [2020 05  1];    % 
       m4 = [2020 06  1];    % 
       m5 = [2020 07  1];    % 
       m6 = [2020 08  1];    % 
       m7 = [2020 09  1];    % 
       m8 = [2020 10  1];    % 
   
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
      M = 'AMJJASO';
      for c = 1:6
        text(5+tm(c+1), ys*yMax,M(c))
      end

   end
