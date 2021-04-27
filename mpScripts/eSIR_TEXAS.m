%eSIR_TEXAS.m

clear
close all
clc

global a b

% Initialize matrices
  num = 5000;
  S  = zeros(num,1);       % Susceptible population
  I = zeros(num,1);        % Active infected population
  R = zeros(num,1);        % Removals form infected population  
%  C = zeros(num,1);       % ReCovered population
%  D = zeros(num,1);       % Dead population
  
% Model Parameters  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  cn = 'TEXAS';
  a = 0.130;
  b = 0.050;
  f = 77e4;
  
  k0 = 1.10e-5;
  D0 = 5e3;
  
  k1 = 0.10e-5; D1 = 28000; ds1 = 1300;
  
  k2 = 0.20e-5; D2 = 15000; ds2 = 2500;
  
  tMax = 500;
  
  I(1) = 1.5e-5;
  R(1) = 0;
  
% Starting Date
  zSTART = datenum([2020 3 12]);  
  
  
% Setup  ==============================================================
  S(1) = 1;
  t = linspace(0,tMax,num);
  h = t(2) - t(1);
  
  
% SOLVE ODEs with RK4  ================================================
for n = 1 : num-1
        
   kS1 = SDOT(t(n), S(n), I(n), R(n));
   kI1 = IDOT(t(n), S(n), I(n), R(n));
   kR1 = RDOT(t(n), S(n), I(n), R(n));
   
   kS2 = SDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   kI2 = IDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   kR2 = RDOT(t(n) + h/2, S(n) + kS1/2, I(n) + kI1/2, R(n) + kR1/2);
   
   kS3 = SDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   kI3 = IDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   kR3 = RDOT(t(n) + h/2, S(n) + kS2/2, I(n) + kI2/2, R(n) + kR2/2);
   
   kS4 = SDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   kI4 = IDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   kR4 = RDOT(t(n) + h,   S(n) + kS3,   I(n) + kI3,   R(n) + kR3);
   
   S(n+1) = S(n) + h*(kS1+2*kS2+2*kS3+kS4)/6;
   I(n+1) = I(n) + h*(kI1+2*kI2+2*kI3+kI4)/6;
   R(n+1) = R(n) + h*(kR1+2*kR2+2*kR3+kR4)/6;
   
    if n > 1880 && n < 2401; S(n+1) = 0.50; end
    if n > 2380 && n < 2651; S(n+1) = 0.56; end
    if n > 2650 && n < 3001; S(n+1) = 0.46; end
    if n > 3000 && n < 3021; S(n+1) = 0.60; end
    
end
  
  I = f.*I;
  R = f.*R;

  
  D = D0.*(1 - exp(-k0.*R));
   
  D(ds1:ds2) = D(ds1) + D1.*(1 - exp(-k1.*(R(ds1:ds2)-R(ds1))));
  D(ds2:end) = D(ds2) + D2.*(1 - exp(-k2.*(R(ds2:end)-R(ds2))));
  
%   D = D0.*(1 - exp(-k0.*R));
%   
%   
%   D(nS:end) = +D(nS) + D1.*(1 - exp(-k1.*(R(nS:end)-R(nS))));
  
  C = R - D;
  Itot = I + R;
  
  Idot = a.*S.*I - b.*I;
  
% Percentages: model
  pI = 100*I(end) / Itot(end);
  pC = 100*C(end) / Itot(end);
  pD = 100*D(end) / Itot(end);
  
% effective reproductive rate
  Re = (a/b).*S;

  
% =====================================================================
% DATA: Id Rd Dd / Number of data days for COVID19 
  covid = sir;
 % Ndays = length(covid);
  Idtot = covid(:,1);
  Id    = covid(:,2);
  Dd    = covid(:,3);
    
  Cd = Idtot - Id - Dd;
  Rd = Cd + Dd;
  Idmax = max(Id);
  Cdmax = max(Cd);
  Ddmax = max(Dd);
  ITdmax = max(Idtot);
  
  Imax = max(I);
  Cmax = max(C);
  Dmax = max(D);
  Smax = max(S);
  ITmax = max(Itot);
  
  Ndays = length(Idtot);
  
% Percentages: Data
  pId = 100*Id(end) / Idtot(end);
  pCd = 100*Cd(end) / Idtot(end);
  pDd = 100*Dd(end) / Idtot(end);

  
EI = 0; Edead = 0; ER = 0; EItot = 0;
  for c = 1 : Ndays
   z = find(t>c,1);
   EItot = EItot + (Itot(z) - Idtot(c))^2;
   EI = EI + (I(z) - Id(c))^2;
   ER = ER + (R(z) - Rd(c))^2;
   Edead = Edead + (D(z) - Dd(c))^2;
  end
  E = sqrt(EItot + EI + ER);
  fprintf('EItot    = %2.2e \n',sqrt(EItot))
  fprintf('EI       = %2.2e \n',sqrt(EI))
  fprintf('ER       = %2.2e \n',sqrt(ER))
  fprintf('Edead    = %2.2e \n',sqrt(Edead))
  fprintf('E        = %2.2e \n',sqrt(E))
  
  
% GRAPHICS  ===========================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.00 0.01 0.45 0.85]);
   set(gcf,'color','w');
   FS = 14;
   
subplot(5,1,1)   
   H = plot(1:Ndays,Idtot,'b+','linewidth',1);
%   set(H,'markersize',3,'markerfacecolor','b')
   
  % xlabel('Days elapsed')
   ylabel('Total  infections')
   title(cn)
   hold on
   plot(t,Itot,'r','linewidth',2)
   xlim([0 tMax])
   set(gca,'xtick',0:50:tMax)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   hold on
   grid on; box on;
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   grid on
 
 subplot(5,1,2) 
   plot(t,I,'r','linewidth',2)
   grid on; box on;
 %  xlabel('Days elapsed')
   ylabel('Active  infections')
   hold on
   H = plot(1:Ndays,Id,'b+','linewidth',1);
 %  set(H,'markersize',3,'markerfacecolor','b')
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   xlim([0 tMax])
   
      
  subplot(5,1,3)   
   H = plot(1:Ndays,Rd,'b+','linewidth',1);
 % set(H,'markersize',3,'markerfacecolor','b')
   grid on; box on;
 %  xlabel('Days elapsed')
   ylabel('Removed')
   hold on
   plot(t,R,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   xlim([0 tMax])
   
   subplot(5,1,4)   
   H = plot(1:Ndays, Cd,'b+','linewidth',1);
 % set(H,'markersize',3,'markerfacecolor','b')
   grid on; box on;
 %  xlabel('Days elapsed')
   ylabel('ReCovered')
   hold on
   plot(t,C,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   xlim([0 tMax])
   
  subplot(5,1,5)   
   H = plot(1:Ndays,Dd,'b+','linewidth',1);
%   set(H,'markersize',3,'markerfacecolor','b')
   grid on; box on;
   xlabel('Days elapsed')
   ylabel('Deaths')
   hold on
   plot(t,D,'r','linewidth',2)
   text(-30,-7e3,'12 Mar 2020')
 %  ytickformat('%2.1e')
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   xlim([0 tMax])
   
   
figure(2)  % =========================================================
   set(gcf,'units','normalized');
   set(gcf,'position',[0.52 0.01 0.40 0.85]);
   set(gcf,'color','w');
   FS = 14;
   
subplot(5,1,1)   
   plot(t,S,'r','linewidth',2)
   xlabel('Days elapsed')
   ylabel('Susceptible  S')
   title(cn)
   xlim([0 tMax])
   ylim([0 1.1])
 %  set(gca,'ytick',0:0.2:1)
   set(gca,'xtick',0:50:tMax)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   hold on
   grid on; box on
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   
   plot([0 tMax],[b/a b/a],'m','linewidth',1)
   txt = sprintf('S_C = %2.2f',b/a);
   Htext = text(400,0.40,txt);
   set(Htext,'fontName','times','fontSize',12,'BackgroundColor','w')
   
   z = find(t >Ndays,1);
   Htext = plot(t(z),S(z),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
   text(-30,-0.3,'12 Mar 2020')
   
subplot(5,1,2)
   plot(t,Idot,'m')
   hold on
   plot(t,gradient(I,h),'r','linewidth',2)
   ylabel('dI/dt')
   grid on
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
% subplot(4,1,1)   
%    plot(t,Re,'r','linewidth',2)
%    xlabel('Days elapsed')
%    ylabel('Eff.  Reproductive  rate')
%    % ylabel('Susceptible  S')
%    title(cn)
%    xlim([0 tMax])
%  %  ylim([0 1.1])
%  %  set(gca,'ytick',0:0.2:1)
%    set(gca,'xtick',0:50:tMax)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    hold on
%    grid on; box on
%    yMax = ylim;
%    z = yMax(2);
%    plotMonths(z)
%    ylim([0 z])
%    plot([0 tMax],[1 1],'m','linewidth',1)
% %    txt = sprintf('S_C = %2.2f',b/a);
% %    Htext = text(150,0.8,txt);
% %    set(Htext,'fontName','times','fontSize',12,'BackgroundColor','w')
%    
%    z = find(t >Ndays,1);
%    Htext = plot(t(z),Re(z),'ro');
%    set(Htext,'markersize',6,'MarkerFaceColor','r')
%    text(-30,-3,'22  Jan 2020')

 subplot(5,1,3)   
   H = plot(Rd,Dd,'b+','linewidth',1);
 %  set(H,'markersize',3,'markerfacecolor','b')
   grid on; box on;
   xlabel('Removals')
   ylabel('Deaths')
   hold on
   plot(R,D,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   
subplot(5,1,4)   
   H = plot(Id,Rd,'b+','linewidth',1);
 %  set(H,'markersize',3,'markerfacecolor','b')
   grid on; box on;
   xlabel('Active infections')
   ylabel('Removals')
   hold on
   plot(I,R,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   
  subplot(5,1,5)   
   H = plot(Idtot,Rd,'b+','linewidth',1);
  % set(H,'markersize',3,'markerfacecolor','b')
   grid on; box on;
   xlabel('Total infections')
   ylabel('Removals')
   hold on
   plot(Itot,R,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   

figure(9)  % ===========================================================
   set(gcf,'units','normalized');
   set(gcf,'position',[0.67 0.10 0.30 0.4]);
   set(gcf,'color','w');
   FS = 14;
   
   xlim([0 120])
   ylim([0 250])
   hh = 250; dh = -20;
   
   txt = cn;
   Htext = text(0,hh,txt,'fontsize',12);
   set(Htext,'color','k')
   txt = 'Start Date';
   Htext = text(30,hh,txt,'fontsize',12);
   set(Htext,'color','k')
   z = zSTART;
   z = datetime(z,'ConvertFrom','datenum');
   z.Format = 'dd-MMM-yyyy';
   txt = cellstr(z) ; 
   text(55,hh,txt,'fontsize',12)
   
   hh = hh+2*dh;
   txt = 'MODEL PARAMETERS';
   text(0,hh,txt,'fontsize',12)
   
    hh = hh+dh;
   txt = sprintf('I(0) = %3.2e  R(0) = %2.2e  f = %2.2e  a = %2.3f  b = %2.3f', I(1)/f, R(1)/f, f,a ,b);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
 % txt = 'D = D_0 (1-exp(-k_0 R))';
 %  text(0,hh,txt,'fontsize',12)
   txt = sprintf('D_0 = %3.2e   k_0 = %3.2e', D0, k0);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+2*dh; 
   txt = 'DATA';
   text(0,hh,txt,'fontsize',12)
   
   z = zSTART;
   z = z + Ndays+3;
   z = datetime(z,'ConvertFrom','datenum');
   z.Format = 'dd-MMM-yyyy';
   txt = cellstr(z) ; 
   text(20,hh,txt,'fontsize',12)
    
   hh = hh+dh;
   z1 = Id(end); z2 = Cd(end); z3 = Dd(end); z4 = z1+z2+z3;  
   txt = sprintf('I = %5.0f    C = %5.0f    D = %5.0f    I_{tot} = %5.0f', z1,z2,z3,z4);
   text(6,hh,txt,'fontsize',12)
   
     hh = hh+dh;
     txt = sprintf('percent: Active = %2.1f  ReCoveries = %2.1f    Deaths = %2.1f   ',pId,pCd,pDd);
     text(2,hh,txt,'fontsize',12)  
   
     hh = hh+2*dh; 
     txt = 'MODEL PREDICTIONS';
     text(0,hh,txt,'fontsize',12)
   
     z = zSTART;
     z = z + tMax;
     z = datetime(z,'ConvertFrom','datenum');
     z.Format = 'dd-MMM-yyyy';
     txt = cellstr(z) ; 
     text(58,hh,txt,'fontsize',12)
   
     hh = hh+dh;
     z1 = I(end); z2 = C(end); z3 = D(end); z4 = z1+z2+z3;  
     txt = sprintf('I = %5.0f    C = %5.0f    D = %5.0f    I_{tot} = %5.0f', z1,z2,z3,z4);
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
     text(60,hh,txt,'fontsize',12)
   
     hh = hh+dh;
     txt = sprintf('percent: Active = %2.1f  ReCoveries = %2.1f    Deaths = %2.1f   ',pI,pC,pD);
     text(2,hh,txt,'fontsize',12)   
        
     axis off
    
%  figure(5)
%    plot(t,Re,'r','linewidth',2)
%    xlabel('Days elapsed')
%    ylabel('Eff.  Reproductive  rate')
%    title(cn)
%    xlim([0 tMax])
%    ylim([0 1.1*max(Re)])
%  %  set(gca,'ytick',0:0.2:1)
%    set(gca,'xtick',0:50:tMax)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    hold on
%    grid on; box on
%    yMax = ylim;
%    z = yMax(2);
%    plotMonths(z)
%    
%    plot([0 tMax],[1  1],'m','linewidth',1)
%  %  txt = sprintf('S_C = %2.2f',b/a);
%  %  Htext = text(210,0.35,txt);
%  %  set(Htext,'fontName','times','fontSize',12,'BackgroundColor','w')
%    
%    z = find(t >Ndays,1);
%    Htext = plot(t(z),Re(z),'ro');
%    set(Htext,'markersize',7,'MarkerFaceColor','r')


% FUNCTIONS  ===========================================================

function FS = SDOT(t,S,I,R)
  global a
    FS = -a*S*I;
end


function FI = IDOT(t,S,I,R)
  global a b
    FI = a*S*I - b*I;
end

function FR = RDOT(t,S,I,R)
  global b
    FR = b*I;
end

% ---------------------------------------------------------------------

function  plotMonths(z)
     
       yMax = z;
     % Months
     
       m1 = [2020 03 12];    % day zero
       m2 = [2020 04 1];      % Feb
       m3 = [2020 05  1];    % Mar
       m4 = [2020 06  1];    % Apr
       m5 = [2020 07  1];    % May
       m6 = [2020 08  1];    % Jun
       m7 = [2020 09  1];    % Jul
       m8 = [2020 10  1];    % Aug
       m9 = [2020 11  1];    % Sep
       m10 = [2020 12  1];    % Oct
       m11 = [2021 1  1];    % Nov
       m12 = [2021 2  1];   % Dec
       m13 = [2021 3  1];   % Jan
       m14 = [2021 4  1];
       m15 = [2021 5  1];
       m16 = [2021 6  1];
       m17 = [2021 7  1];
       m18 = [2021 8  1];
       m19 = [2021 9  1];
       m20 = [2021 10  1];
       m21 = [2021 11  1];   
          
       
   
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
      tm(21) = datenum(m21) - tm(1);
      
      
      for n = 2:18
        plot([tm(n), tm(n)], [0, yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
        hold on
      end
      
      ys = 0.9;
      text(30, ys*yMax,'A')
      text(60, ys*yMax,'M')
      text(90, ys*yMax,'J')
      text(120,ys*yMax,'J')
      text(150,ys*yMax,'A')
      text(180,ys*yMax,'S')
      text(210,ys*yMax,'0')
      text(240,ys*yMax,'N')
      text(270,ys*yMax,'D')
      text(300,ys*yMax,'J')
      text(330,ys*yMax,'F')
      text(360,ys*yMax,'M')
      text(390,ys*yMax,'A')
      text(420,ys*yMax,'M')
      text(450,ys*yMax,'J')
      text(480,ys*yMax,'J')
      
   end





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
T = [314898 155158 3801];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [327532 153970 3981];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [330501 156913 4007];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [339210 162211 4063];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [352165 170059 4235];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [357127 166299 4550];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [375291 166755 4710];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [377396 168808 4762];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [393653 167076 5067];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [393683 167104 5069];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [399137 164853 5177];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [409625 159292 5884];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [414877 164424 6004];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [430252 163240 6470];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [433270 166031 6703];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [448652 158782 7266];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [454364 164419 7341];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [459021 154218 7381 ];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [470195 156306 7627];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [473206 159234 7710];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [488499 156608 8087];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [497406 157394 8344];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [503891 157051 8497];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [508124 154699 8580];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [516774 158265 8676];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [526680 159364 9004];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [527441 160108 9021];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [539989 154742 9487];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [548911 155458 9736];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [556591 153447 9878];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [559810 150236 10003];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [571129 155234 10078];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [578692 152545 10244];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [585454 147652 10517];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [589246 144710 10656];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [596391 139167 11553];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [599386 136075 11717];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [601303 137500 11722];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [607547 133464 11835];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [610288 131267 11878];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [617443 128850 12087];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [625294 117961 12499];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [631118 117736 12667];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [635641 110825 12785];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [636317 111467 12788];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [641950 108729 12882];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [646268 103866 13021];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [647898 103894 13061];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [655057 100629 13414];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [661505 93134 13561];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [665113 88391 13693];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [667650 88266 13777];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [668128 86608 13793];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [671188 87150 13853 ];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [677890 85285 14034];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [680563 84944 14056];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [686477 85501 14345];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [687029 85368 14345];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [691553 83146 14517];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [696697 84039 14536];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [701009 82312 14699];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [704813 83116 14730];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [712946 82489 14951];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [713830 82874 14958];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [719262 83027 15083];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [720281 82670 15172];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [731989 88720 15216];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [733173 88933 15254];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [751571 88780 15361];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [758286 88540 15592];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [760941 87946 15629];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [768203 89146 15823];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [768869 89635 15826];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [774504 84027 15868];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [777165 85489 15923];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [783468 87923 16023];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [787271 87470 16132];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [790467 87069 16251];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [800336 89237 16409];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [802180 90016 16480];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [805485 90373 16505];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [812530 91546 16638];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [814323 91449 16664];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [822424 93398 16857];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [824812 95141 16886];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;% T = [];
T = [831535 96823 17052];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [832928 97758 17055];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [837509 99482 17084];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [839569 100420 17098];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [845656 103073 17196];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [854748 104964 17384];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [857102 106391 17393];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [866811 109844 17552];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [867394 110118 17554];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [872863 113769 17582];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [880867 114527 17664];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [889513 115605 17786];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [891907 115792 17850];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [894607 118100 17878];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [907172 121230 18021];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [909749 123102 18026];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [917873 126457 18093];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [925356 128405 18188];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [934813 132545 18289];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [940761 134409 18392];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [942486 135017 18402];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [955160 138154 18578];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [969247 144589 18662];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [985582 151964 18887];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [988325 153768 18904];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1007155 162469 19154];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1012428 163723 19247];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1019926 162636 19286];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1027114 168800 19324];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1030227 171061 19330];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1055199 180123 19573];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1059222 184474 19691];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1078194 181080 19946];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1085490 177567 20059];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1086987 178292 20064];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1094356 181750 20140];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1106235 186438 20218];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1130292 195493 20540];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1147276 205036 20766];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1161219 210321 20956];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1169076 213647 21085];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1174105 215516 21164];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1178718 218163 21165];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1206567 230617 21390];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1212127 233712 21471];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1235838 237196 21834];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1238759 237280 21915];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1252214 242067 21972];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1265779 242975 22052];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1280720 244672 22294];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [1299915 245995 22505];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [1317457 252523 22729];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1331516 256966 22977];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1341726 257783 23181];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1352230 259727 23252];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1366709 264520 23335];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1381360 264005 23534];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1396552 267199 23821];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1410035 268788 24066];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1423832 249054 24298];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1470079 264681 24489];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1487369 264665 24615];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1500111 263987 24670];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1520038 261978 24890];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1540991 268001 25146];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1556006 269633 25401];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1556006 269633 25401];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1572594 277869 25663];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1586901 277952 25910];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1597473 281004 26019];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1614482 275225 26161];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1635987 278292 26461];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1663327 284550 26957];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1666099 282329 27109];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1671790 281192 27129];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1678480 282938 27186];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1696714 298254 27187];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1715594 289546 27313];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1739416 294004 27584];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1770886 298926 28213];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1782174 301397 28499];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1796622 305042 28594];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1819120 318906 28654];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1830651 320004 28738];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1868896 323248 29028];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1894162 328696 29357];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1906225 331942 29695];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1941212 348638 30190];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1959400 354663 30527];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1975989 364277 30684];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1997795 363939 30795];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2024742 367811 31100];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2050738 374813 31481];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;    
T = [2071261 376930 31863];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2102400 390626 32292];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2111967 388526 32626];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2133490 400284 32839];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2137109 393796 32880];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2155879 378982 33234];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2192450 382355 33702];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2214767 381088 34165];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2233304 384050 34600];
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









end



