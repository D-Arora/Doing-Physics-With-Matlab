% eSIR.m

clear
close all
clc


global a b 

% Initialize matrices
  num = 5000;
  S  = zeros(num,1);       % Susceptible population factor
  I = zeros(num,1);        % Active infected population
  R = zeros(num,1);        % Removals form infected population 
  D = zeros(num,1);
  A = zeros(num,1);
  B = zeros(num,1);
  P = zeros(num,1);
  Q = zeros(num,1);
  
 % Starting Date  z
   zSTART = datenum([2020 04 01]); 
% LOAD DATA  ==========================================================

  load('covid.mat')

% SETUP
  nT = 5000; tMax = 500;
  t = linspace(0,tMax,nT);
  h = t(2)-t(1);

% SELECT COUNTRY  =====================================================
%   disp('  ')
%   disp('SELECT COUNTRY')
%   disp('1 World   2 USA      3 Texas        4 India    5 Italy   6 Germany')
%   disp('7 Iran    8 Israel   9 Australia   10 Korea   11 Japan')
%   disp('  ')
%flagC = input('  Country = ');

flagC = 4;


% COUNTRIES  ==========================================================   
  switch  flagC

 
 case 1  %   WORLD ****************************************************
         cn = "WORLD    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.4;%0.19
         b = 0.02;%0.034;
         f = 19e5;
           
         S(1) = 0.5;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 1200;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 2500;
      
         

case 2   %   U.S.A. **************************************************
         cn = "U.S.A.    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.20;   
         b = 0.015;  
         f = 12.5e5;
         
         S(1) = 0.8;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 1200;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 2500;
         
     
case 3   %   TEXAS ***************************************************
         cn = "TEXAS    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.130;
         b = 0.040;
         f = 65e4;
         
         S(1) = 0.5;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         
         D00 = 0;
         
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 1200;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 2500;

 
 case 4  %   INDIA *********************************************************
         cn = "INDIA    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         f = 12e6;
         S(1) = 0.55;
       %  I(1) = Id(1)/f;
       %  R(1) = Rd(1)/f;
         I(1) = 1.5*1e-3;
         R(1) = 2*10e-6;
         
        
         mD(1) = 1;
         nD(1) = 0;
         
         k(1) = 0.035e-5;
         mD(2) = 1700;
         nD(2) = 100e3;
         
         k(2) = 0.03e-5;
         mD(3) = 3500;
         nD(3) = 90e3;
         k(3) = 0.06e-5;
         
         k(3) = 0.03e-5;
         mD(4) = 5000;
         nD(4) = 50e3;
         k(4) = 0.03e-5;
         
         
         segments = length(k)-1;
        
         a = 0.20;
         b = 0.08;  
         A = a.*ones(num,1);
         B = b.*ones(num,1);
         
         n1 = 800;  n2 = 2001; P(n1:n2) = 3.2e-4;
         n1 = 2000; n2 = 3201; P(n1:n2) = 2.0e-4;
         n1 = 3200; n2 = 4001; P(n1:n2) = 11.0e-4;
        
         m1 = 1800; m2 = 3001; Q(m1:m2) = 3.2e-5;
         
      
    
         
   case 5  %   ITALY ***********************************************
         cn = "ITALY    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.18;
         b = 0.036;
         f = 1.75e5;
         
         S(1) = 0.6;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0.8e4;
         k0 = 1e-5; D0 = 31e3;
         k1 = 0.04e-5; D1 = 105e3; ds1 = 2100;
         k2 = 0.14e-5; D2 = 0; ds2 = 5000;
         
         
   
  case 6 %     GERMANY *****************************************************
         cn = "GERMANY    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.298;
         b = 0.061;
         f = 17e4;
         
         S(1) = 0.5;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         
         D00 = Dd(1);
         
         k0 = 0.6e-5; D0 = 13e3;
         
         k1 = 0.08e-5; D1 = 40e3; ds1 = 2350;
         
         k2 = 0.17e-5; D2 = 35e3; ds2 = 2700;
         
         
  
 case 7  %     IRAN *****************************************************
         cn = "IRAN    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.15;
         b = 0.061;    %0.0397;
         f = 15e4;
         
         S(1) = 0.5;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 5000;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 5000;
         
           
  case 8 % ISRAEL **************************************************
         cn = "ISRAEL    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.25;
         b = 0.066;
         f = 3e4;
         
         S(1) = 0.5;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         k0 =0.5e-5; D0 = 9e3;
         k1 = 0.10e-5; D1 = 28000; ds1 = 5000;
         k2 = 0.15e-5; D2 = 16000; ds2 = 5000;
         
         
       
  case 9  % AUSTRALIA  ====================================================
         cn = "AUSTRALIA    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a =  0.21 ; %0.20560;
         b = 0.035  ;%0.0346;
         f = 4.8e3;
  
         
         D00 = 0;
         k0 = 2.8e-4;  D0 = 120;
         k1 = 0.5e-4;  D1 = 900;  ds1 = 1100;
         k2 = 1.1e-4;  D2 = 650;  ds2 = 1450;
         
                 
         S(1) = 0.25;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
         
         
         
 %   SOUTH KOREA *************************************************************
    case 10
         cn = "SOUTH KOREA    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
          a = 0.42;
          b = 0.0405;
          f = 1.8e4;
         
         S(1) = 0.03;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         k0 = 0.9e-4;  D0 = 425;  
         k1 = 0.5e-4;  D1 = 520;   ds1 = 1600;
         k2 = 0.30e-4;  D2 = 1225;  ds2 = 2600; 
         
         
  
  % JAPAN  =============================================================       
  case 11
         cn = "JAPAN    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid(:,3*flagC - 2) == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         f = 5.2e4;
 
         S(1) = 0.5;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;       
         
         a = 0.38;
         b = 0.07;  
         A = a.*ones(num,1);
         B = b.*ones(num,1);
     
    % Deaths
         mD(1) = 1;
         nD(1) = 0;
         
         k(1) = 2e-5;
         mD(2) = 2500;
         nD(2) = 3e3;
         
          k(2) = 1e-5;
          mD(3) = 3200;
          nD(3) = 4e3;
          
           k(3) = 1e-5;
           mD(4) = 5000;
           nD(4) = 5e3;
           
           k(4) = 1e5;
                  
         segments = length(k)-1;
        
   % Susceptible   
   
        n1 = 700;  n2 =1401; P(n1:n2) = 12e-4;
        n1 = 1400; n2 = 2001; P(n1:n2) = 6e-4;
        n1 = 2000; n2 = 2401; P(n1:n2) = 20e-4;
        n1 = 2400; n2 = 2701; P(n1:n2) = 40e-4;
        n1 = 2700; n2 = 2801; P(n1:n2) = 80e-4;
        n1 = 2800; n2 = 3001; P(n1:n2) = 120e-4;
        n1 = 3000; n2 = 3401; P(n1:n2) = 30e-4;
        n1 = 3400; n2 = 4000; P(n1:n2) = 50e-4;
       
      %  m1 = 2000; m2 = 2801; Q(m1:m2) = 50e-5;
      %  m1 = 2900; m2 = 4001; Q(m1:m2) = 550e-5;
        
      %  A(2950:end) = 0.27;
      %  B(2950:end) = 0.05;
        
  end
 
     
% SOLVE ODE  =========================================================
for n = 1 : num-1
    
 %  [p1(n), p2(n), p3(n), p4(n)] = CV(flagC,n, q1(n),q2(n),q3(n),q4(n)) ; 
   a = A(n);
   b = B(n); 
  
   S(n) = S(n) + P(n);
   I(n) = I(n) + Q(n);
       
   kS1 = EDOT(t(n), S(n), I(n), R(n));
   kI1 = IDOT(t(n), S(n), I(n), R(n));
   kR1 = RDOT(t(n), S(n), I(n), R(n));
   
   kS2 = EDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   kI2 = IDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   kR2 = RDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   
   kS3 = EDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   kI3 = IDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   kR3 = RDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   
   kS4 = EDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   kI4 = IDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   kR4 = RDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   
   S(n+1) = S(n) + h*(kS1+2*kS2+2*kS3+kS4)/6;
   I(n+1) = I(n) + h*(kI1+2*kI2+2*kI3+kI4)/6;
   R(n+1) = R(n) + h*(kR1+2*kR2+2*kR3+kR4)/6;
   
  
   
end


%   Estimate a and b from data ==========================================
%     Change tN and dt for different estimates
%       and view answers in Command Window
    tN =291; dt = 10;
    
    Rdot = ( Rd(tN+dt) - Rd(tN-dt) )/(2*dt);
    bD = Rdot/Id(tN);
    
 % At a peak  a = b/S
    tN = 291;
    aD = b/S(10*tN);
    
    
    

% CALCULATIONS  ========================================================
  
   Sdays = S(1:10:end)';
  
% Population scaling   
   I = f.*I;
   R = f.*R;

% Deaths
  for c = 1 : segments
      D(mD(c): mD(c+1)) = D(mD(c)) + nD(c+1).*(1 - exp(-k(c).*(R(mD(c) : mD(c+1)) - R(mD(c)))));
      
  end
   
   C = R - D;
   Itot = I + R;
 
% Errors  ===============================================================  
  EI = 0; Edead = 0; ER = 0; EItot = 0;
  for c = 1 : Ndays
   z = find(t>c,1);
   EItot = EItot + (Itot(z) - Idtot(c))^2;
   EI = EI + (I(z) - Id(c))^2;
   ER = ER + (R(z) - Rd(c))^2;
   Edead = Edead + (D(z) - Dd(c))^2;
  end
   Emax = sqrt(EItot + EI + ER);
   EItot = sqrt(EItot);
   EI    = sqrt(EI);
   ER    = sqrt(ER);
   Edead    = sqrt(Edead);
  
   zz = zSTART;
   zz = datetime(zz,'ConvertFrom','datenum');
   zz.Format = 'dd-MMM-yyyy';
   fprintf('Start date    %s \n',zz)
   zz = zSTART+ Ndays;
   zz = datetime(zz,'ConvertFrom','datenum');
   zz.Format = 'dd-MMM-yyyy';
 %  fprintf('Current date  %s \n',zz)
 %  disp('   ')
  
%    fprintf('EItot    = %2.2e \n',EItot)
%    fprintf('EI       = %2.2e \n',EI)
%    fprintf('ER       = %2.2e \n',ER)
%    fprintf('Edead    = %2.2e \n',Edead)
%    fprintf('E        = %2.2e \n',Emax)
%  
%    
% GRAPHICS  =========================================================

figure(9)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.2 0.2 0.45 0.55]);
   set(gcf,'color','w');
   FS = 14;

  fig = figure(9); 
  left_color  = [1 0 0];
  right_color = [0 0 1];
  set(fig,'defaultAxesColorOrder',[left_color; right_color]);   
  
  zz = find(t > Ndays,1);

subplot(3,1,1)   
   plot(t,S,'r','linewidth',2)
 
   ylabel('S     ');
   title(cn)
   xlim([0 tMax])
   set(gca,'xtick',0:50:tMax)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   
   yLH = get(gca,'ylabel');
   set(yLH,'rotation',0)
   
   hold on
   grid on; box on
   
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   Htext = plot(t(zz),S(zz),'o');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
       
   
 subplot(3,1,2) 
  %yMax = ylim;
  % z = yMax(2);
   z = 1;
   plotMonths(z)
  hold on
   yyaxis left
   xP = t; yP = A;
   plot(xP,yP,'b','linewidth',2)
   hold on   
   HyL = ylabel('a     ');
   set(HyL,'Color',[0 0 1])
   set(gca,'fontsize',FS,'YColor',[0 0 1])
   yLH = get(gca,'ylabel');
   set(yLH,'rotation',0)
   ylim([0 1])
   grid on; box on 
  
 
   yyaxis right
   xP = t; yP = B;
   plot(xP,yP,'r','linewidth',2)
      
   ylabel('     b')
   set(gca,'fontsize',FS)
   yLH = get(gca,'ylabel');
   set(yLH,'rotation',0)
   ylim([0 1])
   
  
  
 
 subplot(3,1,3) 
   yMax = ylim;
   z = yMax(2);
   z = max(P);
   plotMonths(z)
   hold on
   
   yyaxis left
   xP = t; yP = P;
   plot(xP,yP,'b','linewidth',2)
   hold on   
   ylabel('P     ')
   xlabel('Days elapsed')
   set(gca,'fontsize',FS)
   yLH = get(gca,'ylabel');
   set(yLH,'rotation',0)
   grid on; box on
   
  
   yyaxis right
   xP = t; yP = Q;
   plot(xP,yP,'r','linewidth',2)
      
   ylabel('     Q')
   set(gca,'fontsize',FS)  
   yLH = get(gca,'ylabel');
   set(yLH,'rotation',0)
   
  
   
   
% ---------------------------------------------------------------------
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.00 0.01 0.45 0.85]);
   set(gcf,'color','w');
   FS = 14;
   
subplot(4,1,1)   
   H = plot(1:Ndays,Idtot,'b+','linewidth',1);
   set(H,'markersize',4,'markerfacecolor','b')
   hold on
   plot(t,Itot,'r','linewidth',2)
   yLabel = 'Total  infections';
   graph(tMax, FS, yLabel)
   title(cn)
   
  
 subplot(4,1,2)   
   H = plot(1:Ndays,Id,'b+','linewidth',1);
   set(H,'markersize',4,'markerfacecolor','b')
   hold on
   plot(t,I,'r','linewidth',2)
   yLabel = 'Active  infections';
   graph(tMax, FS, yLabel) 
   
  
  subplot(4,1,3)   
   H = plot(1:Ndays,Cd,'b+','linewidth',1);
   set(H,'markersize',4,'markerfacecolor','b')
   hold on
   plot(t,C,'r','linewidth',2)
   yLabel = '   reCoveries';
   graph(tMax, FS, yLabel) 
   
  subplot(4,1,4)   
   H = plot(1:Ndays,Dd,'b+','linewidth',1);
   set(H,'markersize',4,'markerfacecolor','b')
   hold on
   plot(t,D,'r','linewidth',2)
   yLabel = 'Deaths';
   graph(tMax, FS, yLabel)
   xlabel('Days elapsed')
 
   
   
figure(2)  % =========================================================
   set(gcf,'units','normalized');
   set(gcf,'position',[0.52 0.01 0.40 0.85]);
   set(gcf,'color','w');
   FS = 16;
   
   zz = find(t > Ndays,1);

subplot(4,1,1)   
   plot(t,S,'r','linewidth',2)
 %  area(t,S,'facecolor',[1 0.6 0.6])
   xlabel('Days elapsed')
   ylabel('S');
   title(cn)
   xlim([0 tMax])
 %  ylim([0 1.1])
 %  set(gca,'ytick',0:2:8)
   set(gca,'xtick',0:50:tMax)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   hold on
   grid on; box on
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   
   Htext = plot(t(zz),S(zz),'o');
   set(Htext,'markersize',8,'MarkerFaceColor','r')
     
 subplot(4,1,2)   
   plot(Rd,Dd,'b+','linewidth',1);
   grid on; box on;
   xlabel('Removals')
   ylabel('Deaths')
   hold on
   plot(R,D,'r','linewidth',2)
   
   Htext = plot(Rd(end),Dd(end),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
   
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   
   
subplot(4,1,3)   
   plot(Id,Rd,'b+','linewidth',1);
   grid on; box on;
   xlabel('Active infections')
   ylabel('Removals')
   hold on
   plot(I,R,'r','linewidth',2)
   
   Htext = plot(I(zz),R(zz),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
   
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   

 subplot(4,1,4)   
   H = plot(Idtot,Rd,'b+','linewidth',1);
  % set(H,'markersize',3,'markerfacecolor','b')
   grid on; box on;
   xlabel('Total infections')
   ylabel('Removals')
   hold on
   plot(Itot,R,'r','linewidth',2)
   
   Htext = plot(Itot(zz),R(zz),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
   
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
      
 
   
%  FUNCTIONS  ==========================================================

function graph(tMax, FS, yLabel)
   xlim([0 tMax])
   set(gca,'xtick',0:50:tMax)
   grid on; box on;
   ylabel(yLabel)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
end   


function  plotMonths(z)
       yMax = z;
       m1 = [2020 04  1];    
       m2 = [2020 05  1];    
       m3 = [2020 06  1];    
       m4 = [2020 07  1];   
       m5 = [2020 08  1];    
       m6 = [2020 09  1];    
       m7 = [2020 10  1];   
       m8 = [2020 11  1];  
       m9 = [2020 12  1];   
       m10 = [2021 1  1];   
       m11 = [2021 2  1];
       m12 = [2021 3  1];
       m13 = [2021 4  1];
       m14 = [2021 5  1];
       m15 = [2021 6  1];
       m16 = [2021 7  1];
       m17 = [2021 8  1];
       m18 = [2021 9  1];   
       m19 = [2021 10  1];
       m20 = [2021 11  1];  
          
      tm(1) = datenum(m1);
      tm(2) = datenum(m2) - tm(1);
      tm(3) = datenum(m3) - tm(1);
      tm(4) = datenum(m4) - tm(1);
      tm(5) = datenum(m5) - tm(1);
      tm(6) = datenum(m6) - tm(1);
      tm(7) = datenum(m7) - tm(1);
      tm(8) = datenum(m8) - tm(1);
      tm(9) = datenum(m9) - tm(1);
      tm(10) = datenum(m10) - tm(1);
      tm(11) = datenum(m11) - tm(1);
      tm(12) = datenum(m12) - tm(1);
      tm(13) = datenum(m13) - tm(1);
      tm(14) = datenum(m14) - tm(1);
      tm(15) = datenum(m15) - tm(1);
      tm(16) = datenum(m16) - tm(1);
      tm(17) = datenum(m17) - tm(1);
      tm(18) = datenum(m18) - tm(1);
      tm(19) = datenum(m19) - tm(1);
      tm(20) = datenum(m20) - tm(1);
          
      for n = 2:17
        plot([tm(n), tm(n)], [0, 0.95*yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
        hold on
      end
      
      ys = 0.9;
      
      text(10, ys*yMax,'A')
      text(40, ys*yMax,'M')
      text(70,ys*yMax,'J')
      text(100,ys*yMax,'J')
      text(130,ys*yMax,'A')
      text(160,ys*yMax,'S')
      text(190,ys*yMax,'O')
      text(220,ys*yMax,'N')
      text(250,ys*yMax,'D')
      text(280,ys*yMax,'J')
      text(310,ys*yMax,'F')
      text(340,ys*yMax,'M')
      text(370,ys*yMax,'A')
      text(400,ys*yMax,'M')
      text(430,ys*yMax,'J')
      text(460,ys*yMax,'J')
      text(490,ys*yMax,'A')
end
   
function FS = EDOT(t,S,I,R)
  global a
    FS = -a*S*I;
end


function FI = IDOT(t,S,I,R)
  global a b
    FI =  a*S*I - b*I;
end

function FR = RDOT(t,S,I,R)
  global b
    FR = b*I;
end

