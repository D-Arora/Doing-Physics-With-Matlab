% eSIR.m

clear
close all
clc


global a b

% Initialize matrices
  num = 5000;
  E  = zeros(num,1);       % Effective Reproduction number
  I = zeros(num,1);        % Active infected population
  R = zeros(num,1);        % Removals form infected population 

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

flagC = 7

% COUNTRIES  ==========================================================   
  switch  flagC

 %   WORLD *************************************************************
    case 1
         cn = "WORLD    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.27;%0.19
         b = 0.02;%0.034;
         f = 19e5;
           
         E(1) = 1.3;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 1200;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 2500;
      
         
%   TEXAS *************************************************************
    case 2
         cn = "U.S.A.    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.20;   
         b = 0.015;  
         f = 10.5e5;
         
         E(1) = 1.3;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 1200;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 2500;
         
         
   %   TEXAS *************************************************************
    case 3
         cn = "TEXAS    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.130;
         b = 0.040;
         f = 72e4;
         
         E(1) = 1.3;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         
         D00 = 0;
         
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 1200;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 2500;
         
  
  %   INDIA *************************************************************
    case 4
         cn = "INDIA    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.2;
         b = 0.08;
         f = 12e6;
         
         E(1) = 1.3;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         
         D00 = 0;
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 1200;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 2500; 
         
         
    %   ITALY *************************************************************
    case 5
         cn = "ITALY    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.18;
         b = 0.035;
         f = 1.8e5;
         
         E(1) = 1.8;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0.8e4;
         k0 = 1e-5; D0 = 30e3;
         k1 = 0.04e-5; D1 = 100e3; ds1 = 2100;
         k2 = 0.15e-5; D2 = 0; ds2 = 5000;
         
         
    %   GERMANY *************************************************************
    case 6
         cn = "GERMANY    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.298;
         b = 0.061;
         f = 13.3e4;
         
         E(1) = 1.3;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 1200;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 2500;
         
         
  
    case 7    %  IRAN *****************************************************
         cn = "IRAN    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.15;
         b = 0.061;    %0.0397;
         f = 15e4;
         
         E(1) = 1.3;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         k0 =0.5e-5; D0 = 9e3;
         
         k1 = 0.10e-5; D1 = 28000; ds1 = 5000;
         
         k2 = 0.15e-5; D2 = 16000; ds2 = 5000;
         
           
  case 8     % ISRAEL **************************************************
         cn = "ISRAEL    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.25;
         b = 0.066;
         f = 3e4;
         
         E(1) = 1.8;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         k0 =0.5e-5; D0 = 9e3;
         k1 = 0.10e-5; D1 = 28000; ds1 = 5000;
         k2 = 0.15e-5; D2 = 16000; ds2 = 5000;
         
         
       
  case 9  % AUSTRALIA  ====================================================
         cn = "AUSTRALIA    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a =  0.080 ; %0.20560;
         b = 0.040  ;%0.0346;
         f = 4.8e3;
  
         D00 = 0;
         k0 = 2.8e-4;  D0 = 120;
         k1 = 0.5e-4;  D1 = 900;  ds1 = 1100;
         k2 = 1.1e-4;  D2 = 650;  ds2 = 1450;
         
                 
         E(1) = 0.1;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
         
 %   SOUTH KOREA *************************************************************
    case 10
         cn = "SOUTH KOREA    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
          a = 0.42;
          b = 0.0405;
          f = 1.65e4;
         
         E(1) = 1.3;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;
  
         D00 = 0;
         k0 = 0.9e-4;  D0 = 425;  
         k1 = 0.5e-4;  D1 = 520;   ds1 = 1600;
         k2 = 0.30e-4;  D2 = 1225;  ds2 = 2600;     
  
  % JAPAN  =============================================================       
  case 11
         cn = "JAPAN    1 Apr 20  to  13 Aug 21";
         
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
         a = 0.262;
         b = 0.07;  %0.068;
         f = 3.7e4;
 
         D00 = 0;
         k0 = 4.0e-5;  D0 = 1.3e3;
         k1 = 0.9e-5; D1 = 2.5e3;  ds1 = 1500;
         k2 = 0.5e-5; D2 = 5e3;   ds2 = 2500;
         
         E(1) = 3;
         I(1) = Id(1)/f;
         R(1) = Rd(1)/f;       
              
  end
  
 
 % SOLVE ODE  ========================================================= 
for n = 1 : num-1
        
   kS1 = EDOT(t(n), E(n), I(n), R(n));
   kI1 = IDOT(t(n), E(n), I(n), R(n));
   kR1 = RDOT(t(n), E(n), I(n), R(n));
   
   kS2 = EDOT(t(n) + h/2, E(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   kI2 = IDOT(t(n) + h/2, E(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   kR2 = RDOT(t(n) + h/2, E(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   
   kS3 = EDOT(t(n) + h/2, E(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   kI3 = IDOT(t(n) + h/2, E(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   kR3 = RDOT(t(n) + h/2, E(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   
   kS4 = EDOT(t(n) + h,   E(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   kI4 = IDOT(t(n) + h,   E(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   kR4 = RDOT(t(n) + h,   E(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   
   E(n+1) = E(n) + h*(kS1+2*kS2+2*kS3+kS4)/6;
   I(n+1) = I(n) + h*(kI1+2*kI2+2*kI3+kI4)/6;
   R(n+1) = R(n) + h*(kR1+2*kR2+2*kR3+kR4)/6;
   
   E(n+1) = ERN(flagC,E(n+1),n) ; 
end


% CALCULATIONS  ========================================================
  
   I = f.*I;
   R = f.*R;

   D = D00 + D0.*(1 - exp(-k0.*R));
   D(ds1:ds2) = D(ds1) + D1.*(1 - exp(-k1.*(R(ds1:ds2)-R(ds1))));
   D(ds2:end) = D(ds2) + D2.*(1 - exp(-k2.*(R(ds2:end)-R(ds2))));
 
  C = R - D;
  Itot = I + R;
   
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
   fprintf('Current date  %s \n',zz)
   disp('   ')
  
   fprintf('EItot    = %2.2e \n',EItot)
   fprintf('EI       = %2.2e \n',EI)
   fprintf('ER       = %2.2e \n',ER)
   fprintf('Edead    = %2.2e \n',Edead)
   fprintf('E        = %2.2e \n',Emax)
 
   
  % GRAPHICS  =========================================================
  
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
   H = plot(1:Ndays,Rd,'b+','linewidth',1);
   set(H,'markersize',4,'markerfacecolor','b')
   hold on
   plot(t,R,'r','linewidth',2)
   yLabel = 'R reCoveries';
   graph(tMax, FS, yLabel) 
   
%    %%
%    mx = 8600;
%    xx = 100:150;
%    bx = Rd(133) - mx*133;
%    yy = mx.*xx + bx;
%    
%    plot(xx,yy)
%    
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
   
   zz = find(t >Ndays,1);

subplot(4,1,1)   
   plot(t,E,'r','linewidth',2)
   xlabel('Days elapsed')
   ylabel('S');
   title(cn)
   xlim([0 tMax])
 %  ylim([0 1.1])
   set(gca,'ytick',0:2:8)
   set(gca,'xtick',0:50:tMax)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   hold on
   grid on; box on
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   
   plot([0 tMax],[1 1],'m','linewidth',1)

   Htext = plot(t(zz),E(zz),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
     

 subplot(4,1,2)   
   H = plot(Rd,Dd,'b+','linewidth',1);
 %  set(H,'markersize',3,'markerfacecolor','b')
   grid on; box on;
   xlabel('Removals')
   ylabel('Deaths')
   hold on
   plot(R,D,'r','linewidth',2)
   
   Htext = plot(R(zz),D(zz),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
   
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   
   
subplot(4,1,3)   
   H = plot(Id,Rd,'b+','linewidth',1);
 %  set(H,'markersize',3,'markerfacecolor','b')
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
        plot([tm(n), tm(n)], [0, yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
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
   
function FS = EDOT(t,E,I,R)
  global a
    FS = -a*E*I;
end


function FI = IDOT(t,E,I,R)
  global a b
    FI =  b*(E-1)*I ;%a*E*I - b*I;
end

function FR = RDOT(t,E,I,R)
  global b
    FR = b*I;
end



function z = ERN(flagC,V,n)
      z = V;
switch flagC

case 1   % WORLD ************************************************** 
   if n > 200 && n < 1001;  z = 2.0; end
   if n > 1100 && n < 1301; z = 1.3; end
   if n > 1800 && n < 2401; z = 1.4; end
   if n > 2500 && n < 2901; z = 1.3; end
%  if n > 3000 && n < 3021; z = 1.3; end


case 2  % U.S.A.  *********************************************** 
      if n > 700  && n < 1251;   z = 2.5; end
      if n > 1500  && n < 2000;  z = 0.5; end
%       if n > 1300 && n < 1801;  z = 0.2; end
%       if n > 1800 && n < 2601;  z = 1.0; end
%       if n > 2600 && n < 3001;  z = 1.2; end  


case 3   % TEXAS ************************************************** 
   if n > 200 && n < 1001;  z = 2.0; end
   if n > 1100 && n < 1301; z = 1.3; end
   if n > 1800 && n < 2401; z = 1.4; end
   if n > 2500 && n < 2901; z = 1.3; end
%  if n > 3000 && n < 3021; z = 1.3; end


case 4   % INDIA ************************************************** 
   if n > 200 && n < 1001;  z = 2.0; end
   if n > 1100 && n < 1301; z = 1.3; end
   if n > 1800 && n < 2401; z = 1.4; end
   if n > 2500 && n < 2901; z = 1.3; end
%  if n > 3000 && n < 3021; z = 1.3; end


case 5  % ITALY  *********************************************** 
       if n > 1  && n < 100;   z = 1.8; end
       if n > 1450  && n < 2401;  z = 3.0; end
       if n > 2400 && n < 2701;  z = 0.8; end
       if n > 2700 && n < 2901;  z = 1.5; end
%       if n > 2400 ;  z = 0.8; end


case 6  % GERMANY  *********************************************** 
      if n > 700  && n < 1251;   z = 2.5; end
      if n > 1500  && n < 2000;  z = 0.5; end
%       if n > 1300 && n < 1801;  z = 0.2; end
%       if n > 1800 && n < 2601;  z = 1.0; end
%       if n > 2600 && n < 3001;  z = 1.2; end


case 7   % IRAN ************************************************** 
     if n > 400 && n < 601;  z = 1.5; end
     if n > 800 && n < 1001;  z = 1.3; end
     if n > 1000 && n < 1201;  z = 0.9; end
     if n > 1200 && n < 1401; z = 0.9; end
     if n > 1400 && n < 1601; z = 1.1; end
     if n > 1600 && n < 1801; z = 1.2; end
     if n > 1800 && n < 2001; z = 1.5; end
     if n > 2000 && n < 2201; z = 1.2; end
     if n > 2200 && n < 2401; z = 1.8; end
     if n > 2400 && n < 2601; z = 1.2; end
     if n > 2600 && n < 2801; z = 0.5; end
     if n > 2800 && n < 3001; z = 1.2; end
%      if n > 3000 && n < 3201; z = 1.2; end
%      if n > 3200 && n < 3401; z = 1.2; end
%      if n > 3600 && n < 3801; z = 1.2; end
%      if n > 3800 && n < 4001; z = 1.2; end
%     if n > 4000 && n < 4201; z = 1.2; end
     
     
case 8   % ISRAEL ************************************************** 
    if n > 680 && n < 1151;  z = 2.34; end
    if n > 1300 && n < 1501; z = 0.5; end
    if n > 1500 && n < 1801; z = 2.3; end
    if n > 1800 && n < 2301; z = 0.8; end
 %   if n > 2400 && n < 2901; z = 2.3; end


case 9  % AUSTRALIA  *********************************************** 
      if n > 700  && n < 1251;   z = 2.5; end
      if n > 1500  && n < 2000;  z = 0.5; end
%       if n > 1300 && n < 1801;  z = 0.2; end
%       if n > 1800 && n < 2601;  z = 1.0; end
%       if n > 2600 && n < 3001;  z = 1.2; end


case 10   % SOUTH KOREA ********************************************** 
    if n > 200 && n < 1001;  z = 0.5; end
    if n > 1260 && n < 1401; z = 3.9; end
    if n > 1800 && n < 2101; z = 0.6; end
    if n > 2100 && n < 2761; z = 2.2; end
    if n > 2760; z = 0.7; end


case 11  % JAPAN  *****************************************************
      if n > 900  && n < 1221;   z = 2.5; end
      if n > 1400  && n < 1501;  z = 0.8; end
      if n > 1800 && n < 2001;   z = 1.2; end
      if n > 2100 && n < 2401;   z = 1.8; end
      if n > 2400 && n < 2701;   z = 1.3; end
      if n > 2700 && n < 2951;   z = 1.8; end
      if n > 2950 && n < 3100;   z = 0.8; end
end     
     
end




