%eSIR_INDIA.m

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
  cn = 'INDIA';
  a = 0.2;
  b = 0.08;
  f = 12e6;
  
  k0 = 1.1e-7;
  D0 = 22e4;
  
  tMax = 500;
  
  I(1) = 1e-3;
  R(1) = 10e-6;
  
  
% Starting Date
  zSTART = datenum([2020 3 14]);  
  
  
% Setup  ==============================================================
  S(1) = 0.55;
  t = linspace(0,tMax,num);
  h = t(2) - t(1);
  
  
% SOLVE ODEs with RK4  ================================================
for n = 1 : num-1
    
  if n > 1000 && n < 1300; S(n) = 0.55; end 
  if n > 1600 && n < 1800; S(n) = 0.50; end 
 
  if n > 2000; S(n) = 0.32; end
    
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
   

   
end
  
  I = f.*I;
  R = f.*R;
  
  D = D0.*(1 - exp(-k0.*R));
  
  C = R - D;
  Itot = I + R;
  
  Idot = a.*S.*I - b.*I;
  
% Percentages: model
  pI = 100*I(end) / Itot(end);
  pC = 100*C(end) / Itot(end);
  pD = 100*D(end) / Itot(end);
  
Re = (a/b).*S;




% =====================================================================
% DATA: Id Rd Dd / Number of data days for covid19 
  load('covid19.mat')
  Ndays  = length(covid19(:,1));  
  Idtot  = covid19(53:Ndays,13);
  Id     = covid19(53:Ndays,14);
  Dd     = covid19(53:Ndays,15);  
  
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
  E = (EItot + EI + ER);
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
   text(-30,-45e3,'14 Mar 2020')
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
   Htext = text(50,0.40,txt);
   set(Htext,'fontName','times','fontSize',12,'BackgroundColor','w')
   
   z = find(t >Ndays,1);
   Htext = plot(t(z),S(z),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
   text(-30,-0.6,'14 Mar 2020')
   
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
   z = z + Ndays;
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
     
       m1 = [2020 03 14];    % day zero
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
