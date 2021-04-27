%eSIR_IRAN.m

clear
close all
clc


global a b
% Model Parameters
cn ='Iran';
  a = 0.21;
  b = 0.0397;
  f = 11.5e4;
  
  k0 = 1e-5;
  D0 = 13e3;
  
  k1 = 0.22e-5;
  D1 = 54.5e3;
  ns = 1500;
   
  num = 5000;
  tMax = 500;
  
% Starting Date
  zSTART = datenum([2020 2 26]);  
  
% Initialize matrices
  S  = zeros(num,1);       % Susceptible population
  I = zeros(num,1);        % Active infected population
  R = zeros(num,1);        % Removals form infected population  
%  C = zeros(num,1);        % ReCovered population
%  D = zeros(num,1);        % Dead population
  
  
  S(1) = 1;
  I(1) = 2e-3;
  t = linspace(0,tMax,num);
  h = t(2) - t(1);
  
  
% SOLVE ODEs with RK4  ================================================
for n = 1 : num-1
   
%     if n == 1700
%       a = 0.105;
%       b = 0.076;%0.075;
%       S(n) = 1.1;
%    end 
    
%      if n == 500;  S(n) = 0.2; end
%      if n == 1000; S(n) = 0.8; end 
%      if n == 1250; S(n) = 0.8; end 
%      if n == 1500; S(n) = 0.8; end  
%     if n == 1020; S(n) = 0.5; end  
%     if n == 1080; S(n) = 0.5; end 

     if n > 500 && n < 2900;  S(n) = 0.2; end 
     
     if n > 2250 && n < 2800; S(n) = 0.37; end
%      
%  %   if n > 2400 && n < 2801; S(n) = 0.3; b = 0.03; end
%  %   if n > 2500 && n < 2801; S(n) = 0.60; end
      if n > 2900; S(n) = 0.15;  end
   
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
  
  
  D(ns:end) = D(ns) + D1.*(1 - exp(-k1.*(R(ns:end)-R(ns))));
  
  C = R - D;
  Itot = I + R;
  
  Idot = a.*S.*I - b.*I;
  
  Re = (a/b).*S;

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
  
% Percentages: model
  pI = 100*I(end) / Itot(end);
  pC = 100*C(end) / Itot(end);
  pD = 100*D(end) / Itot(end);
  
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
  EItot = sqrt(EItot);
  EI    = sqrt(EI);
  ER    = sqrt(ER);
  Edead    = sqrt(Edead);
  fprintf('EItot    = %2.2e \n',EItot)
  fprintf('EI       = %2.2e \n',EI)
  fprintf('ER       = %2.2e \n',ER)
  fprintf('Edead    = %2.2e \n',Edead)
  fprintf('E        = %2.2e \n',E)
  
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
   text(-30,-17000,'26 Feb 2020')
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
   text(-30,-0.5,'26 Feb 2020')
   
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
   
%    hh = hh+dh;
%    txt = 'D = D_0 (1-exp(-k_0 R))';
%    text(0,hh,txt,'fontsize',12)
%    txt = sprintf('D_0 = %3.2e   k_0 = %3.2e', D0, k0);
%    text(50,hh,txt,'fontsize',12)
   
   hh = hh+2*dh; 
   txt = 'DATA';
   text(0,hh,txt,'fontsize',12)
   
   z = zSTART;
   z = z + Ndays+2;
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

%      hh = hh+dh;
%      txt = 'Peak';
%      text(10,hh,txt,'fontsize',12)
%      z = find(I == Imax,1);
%      z = t(z) + zSTART;
%      z = datetime(z,'ConvertFrom','datenum');
%      z.Format = 'dd-MMM-yyyy';
%      txt = cellstr(z) ; 
%      text(30,hh,txt,'fontsize',12)
%    
%      txt = sprintf('I_{peak} = %5.0f',Imax);
%      text(60,hh,txt,'fontsize',12)
%    
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
     
       m1 = [2020 02 26];    % day zero
       m2 = [2020 03  1];    
       m3 = [2020 04  1];    
       m4 = [2020 05  1];    
       m5 = [2020 06  1];    
       m6 = [2020 07  1];   
       m7 = [2020 08  1];    
       m8 = [2020 09  1];    
       m9 = [2020 10  1];   
       m10 = [2020 11  1];  
       m11 = [2020 12  1];   
       m12 = [2021 1  1];   
       m13 = [2021 2  1];
       m14 = [2021 3  1];
       m15 = [2021 4  1];
       m16 = [2021 5  1];
       m17 = [2021 6  1];
       m18 = [2021 7  1];
       m19 = [2021 8  1];
       m20 = [2021 9  1];   
       m21 = [2021 10  1];  
          
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
      
      text(24, ys*yMax,'M')
      text(54, ys*yMax,'A')
      text(84, ys*yMax,'M')
      text(114,ys*yMax,'J')
      text(144,ys*yMax,'J')
      text(174,ys*yMax,'A')
      text(204,ys*yMax,'S')
      text(234,ys*yMax,'O')
      text(264,ys*yMax,'N')
      text(294,ys*yMax,'D')
      text(324,ys*yMax,'J')
      text(354,ys*yMax,'F')
      text(384,ys*yMax,'M')
      text(414,ys*yMax,'A')
      text(444,ys*yMax,'M')
      text(474,ys*yMax,'J')
      text(504,ys*yMax,'J')
      
   end


% figure(1)
%    set(gcf,'units','normalized');
%    set(gcf,'position',[0.0 0.10 0.35 0.5]);
%    set(gcf,'color','w');
%    FS = 14;
%    
% subplot(2,2,1)   
%    plot(1:Ndays,Idtot,'b+','linewidth',1)
%    xlabel('Days elapsed')
%    ylabel('Total  infections')
%    title(cn)
%    hold on
%    plot(t,Itot,'r','linewidth',2)
%    xlim([0 tMax])
%    set(gca,'xtick',0:50:tMax)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    hold on
%    grid on; box on;
%    yMax = ylim;
%    z = yMax(2);
%    plotMonths(z)
%    grid on
%  
%  subplot(2,2,2) 
%    plot(1:Ndays,Id,'b+','linewidth',1)
%    grid on; box on;
%    xlabel('Days elapsed')
%    ylabel('Active  infections')
%    hold on
%    plot(t,I,'r','linewidth',2)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    yMax = ylim;
%    z = yMax(2);
%    plotMonths(z)
%    xlim([0 tMax])
%    set(gca,'xtick',0:50:tMax)
%   
%   subplot(2,2,3)   
%    plot(1:Ndays,Cd,'b+','linewidth',1)
%    grid on; box on;
%    xlabel('Days elapsed')
%    ylabel('Recoveries')
%    hold on
%    plot(t,C,'r','linewidth',2)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    yMax = ylim;
%    z = yMax(2);
%    plotMonths(z)
%    xlim([0 tMax])
%    set(gca,'xtick',0:50:tMax)
%    
%   subplot(2,2,4)   
%    plot(1:Ndays,Dd/1000,'b+','linewidth',1)
%    grid on; box on;
%    xlabel('Days elapsed')
%    ylabel('Deaths/1000')
%    hold on
%    plot(t,D/1000,'r','linewidth',2)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    yMax = ylim;
%    z = yMax(2);
%    plotMonths(z)
%     xlim([0 tMax])
%    set(gca,'xtick',0:50:tMax)
%    
% figure(2)  % =========================================================
%    set(gcf,'units','normalized');
%    set(gcf,'position',[0.35 0.10 0.35 0.5]);
%    set(gcf,'color','w');
%    FS = 14;
%    
%    subplot(2,2,1)   
%    plot(t,Re,'r','linewidth',2)
%    xlabel('Days elapsed')
%  %  ylabel('Susceptible  S')
%    ylabel('Reproduction  No.   R_e')
%    title(cn)
%    xlim([0 tMax])
%  %  ylim([0 1.1])
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
%  %  plot([0 tMax],[b/a b/a],'m','linewidth',1)
%    plot([0 tMax],[1 1],'m','linewidth',1)
%   % txt = sprintf('S_C = %2.2f',b/a);
%   % Htext = text(110,1.6*b/a,txt);
%   % set(Htext,'fontName','times','fontSize',12,'BackgroundColor','w')
%    
%    z = find(t >Ndays,1);
%    Htext = plot(t(z),Re(z),'ro');
%    set(Htext,'markersize',6,'MarkerFaceColor','r')
%    
% % subplot(2,2,1)   
% %    plot(t,S,'r','linewidth',2)
% %    xlabel('Days elapsed')
% %    ylabel('Susceptible  S')
% %    title(cn)
% %    xlim([0 tMax])
% %    set(gca,'xtick',0:50:tMax)
% %    ylim([0 1.1])
% %  %  set(gca,'ytick',0:0.2:1)
% %    set(gca,'xtick',0:50:tMax)
% %    set(gca,'fontsize',FS)
% %    set(gca,'fontName','times')
% %    hold on
% %    grid on; box on
% %    yMax = ylim;
% %    z = yMax(2);
% %    plotMonths(z)
% %    
% %    plot([0 tMax],[b/a b/a],'m','linewidth',1)
% %    txt = sprintf('S_C = %2.2f',b/a);
% %    Htext = text(128,0.85,txt);
% %    set(Htext,'fontName','times','fontSize',12,'BackgroundColor','w')
% %    
% %    z = find(t >Ndays,1);
% %    Htext = plot(t(z),S(z),'ro');
% %    set(Htext,'markersize',6,'MarkerFaceColor','r')
% %    
%  subplot(2,2,2)   
%    plot(Rd,Dd,'b+','linewidth',1)
%    grid on; box on;
%    xlabel('Removals')
%    ylabel('Deaths')
%    hold on
%    plot(R,D,'r','linewidth',2)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    
% subplot(2,2,3)   
%    plot(Id,Rd,'b+','linewidth',1)
%    grid on; box on;
%    xlabel('Active infections')
%    ylabel('Removals')
%    hold on
%    plot(I,R,'r','linewidth',2)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    
%   subplot(2,2,4)   
%    plot(Idtot,Rd,'b+','linewidth',1)
%    grid on; box on;
%    xlabel('Total infections')
%    ylabel('Removals')
%    hold on
%    plot(Itot,R,'r','linewidth',2)
%    set(gca,'fontsize',FS)
%    set(gca,'fontName','times')
%    
%    
%    
% figure(9)
%    set(gcf,'units','normalized');
%    set(gcf,'position',[0.67 0.10 0.30 0.4]);
%    set(gcf,'color','w');
%    FS = 14;
%    
%    xlim([0 120])
%    ylim([0 250])
%    hh = 250; dh = -20;
%    
%    txt = 'IRAN:   START DATE   26 February 2020';
%    Htext = text(0,hh,txt,'fontsize',12);
%    set(Htext,'color','k')
%    
%    hh = hh+2*dh;
%    txt = 'MODEL PARAMETERS';
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = sprintf('I(0) = %3.2e  R(0) = %2.2e  f = %2.2e  a = %2.3f  b = %2.3f', I(1)/f, R(1)/f, f,a ,b);
%    text(0,hh,txt,'fontsize',12)
%    
%    hh = hh+dh;
%    txt = 'D = D_0 (1-exp(-k_0 R))';
%    text(0,hh,txt,'fontsize',12)
%    txt = sprintf('D_0 = %3.2e   k_0 = %3.2e', D0, k0);
%    text(50,hh,txt,'fontsize',12)
%    
%    hh = hh+2*dh; 
%    txt = 'DATA';
%    text(0,hh,txt,'fontsize',12)
%    
%    z = zSTART;
%    z = z + Ndays+2;
%    z = datetime(z,'ConvertFrom','datenum');
%    z.Format = 'dd-MMM-yyyy';
%    txt = cellstr(z) ; 
%    text(20,hh,txt,'fontsize',12)
%     
%    hh = hh+dh;
%    z1 = Id(end); z2 = Cd(end); z3 = Dd(end); z4 = z1+z2+z3;  
%    txt = sprintf('I = %5.0f    C = %5.0f    D = %5.0f    I_{tot} = %5.0f', z1,z2,z3,z4);
%    text(6,hh,txt,'fontsize',12)
%    
%    % Recoveries and Deaths
%      hh = hh+dh;
%      txt = sprintf('percent: Active = %2.1f  ReCoveries = %2.1f    Deaths = %2.1f   ',pId,pCd,pDd);
%      text(2,hh,txt,'fontsize',12)  
%    
%      hh = hh+2*dh; 
%      txt = 'MODEL PREDICTIONS';
%      text(0,hh,txt,'fontsize',12)
%    
%      z = zSTART;
%      z = z + tMax;
%      z = datetime(z,'ConvertFrom','datenum');
%      z.Format = 'dd-MMM-yyyy';
%      txt = cellstr(z) ; 
%      text(58,hh,txt,'fontsize',12)
%    
%      hh = hh+dh;
%      z1 = I(end); z2 = C(end); z3 = D(end); z4 = z1+z2+z3;  
%      txt = sprintf('I = %5.0f    C = %5.0f    D = %5.0f    I_{tot} = %5.0f', z1,z2,z3,z4);
%      text(6,hh,txt,'fontsize',12)
% 
%      hh = hh+dh;
%      txt = 'Peak';
%      text(10,hh,txt,'fontsize',12)
%      z = find(I == Imax,1);
%      z = t(z) + zSTART;
%      z = datetime(z,'ConvertFrom','datenum');
%      z.Format = 'dd-MMM-yyyy';
%      txt = cellstr(z) ; 
%      text(30,hh,txt,'fontsize',12)
%    
%      txt = sprintf('I_{peak} = %5.0f',Imax);
%      text(60,hh,txt,'fontsize',12)
%    
%      hh = hh+dh;
%      txt = sprintf('percent: Active = %2.1f  ReCoveries = %2.1f    Deaths = %2.1f   ',pI,pC,pD);
%      text(2,hh,txt,'fontsize',12)   
%         
%      axis off
%     
%      
% 
% 
% % FUNCTIONS  ===========================================================
% 
% function FS = SDOT(t,S,I,R)
%   global a
%     FS = -a*S*I;
% end
% 
% 
% function FI = IDOT(t,S,I,R)
%   global a b
%     FI = a*S*I - b*I;
% end
% 
% function FR = RDOT(t,S,I,R)
%   global b
%     FR = b*I;
% end
% 
% 
% function  plotMonths(z)
%      
%        yMax = z;
%      % Months
%        m1 = [2020 02 26];    % day zero
%        m2 = [2020 03  1];    % 
%        m3 = [2020 04  1];    % 
%        m4 = [2020 05  1];    % 
%        m5 = [2020 06  1];    % 
%        m6 = [2020 07  1];    % 
%        m7 = [2020 08  1];    % 
%        m8 = [2020 09  1];    % 
%        m9 = [2020 10  1];    % 
%        m10 = [2020 11  1];    % 
%        m11 = [2020 12  1];    % 
%        
%        
%    
%       tm(1) = datenum(m1);
%       tm(2) = datenum(m2) - tm(1);
%       tm(3) = datenum(m3) - tm(1);
%       tm(4) = datenum(m4) - tm(1);
%       tm(5) = datenum(m5) - tm(1);
%       tm(6) = datenum(m6) - tm(1);
%       tm(7) = datenum(m7) - tm(1);
%       tm(8) = datenum(m8) - tm(1);
%       tm(9) = datenum(m9) - tm(1);
%       tm(10) = datenum(m10) - tm(1);
%       tm(11) = datenum(m11) - tm(1);
%       
%    
%       for n = 2:11
%         plot([tm(n), tm(n)], [0, yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
%         hold on
%       end
%       
%       ys = 0.94;
%       M = 'MAMJJASOND';
%       for c = 1:length(M)
%         X = 15+(c-1)*30;  
%         text(X, ys*yMax,M(c))
%       end
% 
% end



% FUNCTIONS ============================================================

function covid = sir

c = 1;

T = [139 95 19];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [245 165 26];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [388 281 34];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [593 427 43];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [978 749 54];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [1501 1144 66];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2336 1824 77];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [2922 2278 92];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [3513 2666 108];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [4747 3710 124];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [5823 4009 145];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [6566 4238 194];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [7161 4530 237];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8042 5020 291];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [9000 5687 354];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10075 6370 429];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11364 7321 514];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12729 7779 611];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [13938 8624 724];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14991 9142 853];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [16169 9792 988];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [17361 10516 1135];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [18407 11144 1284];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [19644 11466 1433];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [20610 11419 1556];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [21638 12040 1685];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23049 12861 1812];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [24811 13964 1934];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27017 15315 2077];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [29406 16715 2234];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [32332 18821 2378];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [35408 21212 2517];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [38309 23278 2640];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [41495 24827 2757];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [44605 27051 2898];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [47593 29084 3036];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [50468 30597 3160];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [53183 31954 3294];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [55743 32555 3452];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [58226 32612 3603];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [60500 32525 3739];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [62589 31678 3872];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [64586 30781 3993];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [66220 29801 4110];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [68192 28495 4232];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [70029 23725 4357 ];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [71686 23318 4474];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [73303 22735 4585];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [74877 22065 4683];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [76339 21679 4777];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [77995 20897 4869];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [79494 20472 4958];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [80868 19850 5031];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [82211 20070 5118];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [83505 19023 5209];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [84802 18540 5297];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [85996 17492 5391];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [87026 16702 5481];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [88194 16021 5574];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [89328 15485 5650];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [90481 15114 5710];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [91472 14733 5806];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [92584 14268 5877];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [93657 13909 5957];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [94640 13509 6028];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [95646 13237 6091];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [96448 12942 6156];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [97424 12799 6203];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [98647 12991 6277];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [99970 13155 6340];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [101650 13645 6418];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [103135 13905 6456];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [104691 14313 6541];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [106220 14567 6589];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [107603 14820 6640];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [109286 15179 6685];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [110767 15677 6733];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [112725 16514 6783];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [114533 17140 6854];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [116635 17897 6902];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [118392 18308 6937];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [120198 18746 6988];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [122492 19774 7057];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [124603 20311 7119];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [126949 20958 7183];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [129341 21528 7349];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [131652 22076 7300];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [133521 22090 7359];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [135701 22483 7417];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [137724 22560 7451];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [139511 22566 7508];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [141591 22851 7564];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [143849 23234 7627];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [146668 24060 7677];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [148950 24389 7737];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [151466 24821 7797];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [154445 25563 7878];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [157562 26543 7942 ];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [160696 27478 8012];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [164270 28714 8071];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [167156 29281 8134];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [169425 29178 8209];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [171789 29159 8281];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [173832 29121 8351];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [175927 29045 8425];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [177938 28842 8506];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [180156 28909 8584];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [182525 29217 8659];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [184955 29477 8730];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [187427 29916 8837];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [189876 30336 8950];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [192439 30699 9065];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [195051 31054 9185];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [197647 31384 9272];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [200262 31678 9392];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [202584 31693 9507];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [204952 31738 9623];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [209970 30947 9863];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [212501 30409 9996];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [215096 29863 10130];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [217724 29633 10239];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [220180 29155 10364];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [222669 28851 10508];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [225205 28355 10670];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [227662 28087 10817];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [230211 27766 10958];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [232863 27659 11106];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [235429 27723 11260];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [237878 27521 11408];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [240438 27537 11571];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [243051 27237 11731];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [245688 26757 11931];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [250458 25977 12305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [252720 25258 12447];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [255117 24816 12635];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [257303 24481 12829];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [259652 24081 13032];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [264561 23590 13410];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [267061 22845 13608];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [269440 22776 13791];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [271606 22327 13979];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [273788 21812 14188];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [276202 21710 14405];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [278827 21842 14634];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [281413 21720 14853];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [284034 21730 15074];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [286523 22022 15289];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [288839 22036 15484];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [291172 22259 15700];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [293606 22550 15912];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [296273 23107 16147];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [298909 23450 16343];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [301530 23761 16569];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [304204 23919 16766];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [306752 23940 16982];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [309437 24165 17190];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [312035 24405 17405];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [314786 24634 17617];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [317483 24749 17802];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [320117 24678 17976];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [322567 24711 18132];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [324692 24306 18264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [326712 23914 18427];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [328844 23586 18616];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [331189 23769 18800];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [333699 24467 18988];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [336324 25104 19162];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [338825 25683 19331];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [341070 25948 19492];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [343203 26078 19639];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [345450 26486 19804];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [347835 26982 19972];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [350279 27626 20125];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [352558 28058 20264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [354764 28522 20376];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [356792 28588 20502];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [358906 28798 20643];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [361150 29009 20776];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [363363 29404 20901];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [365606 29716 21020];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [367796 30021 21137];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [369911 30392 21249];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [371816 30610 21359];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [373570 30687 21462];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [375212 30408 21571];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [376894 30098 21672];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [378752 30154 21797];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [380746 30225 21926];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [382772 30420 22044];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [384666 30381 22154];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [386658 30465 22293];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [388810 30828 22410];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [391112 31156 22542];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [393425 31645 22669];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [395488 31848 22798];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [397801 32349 22913];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [399940 32395 23029];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [402029 32630 23157];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [404648 33322 23313];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [407353 33916 23453];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [410334 34683 23632];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [413149 35493 23808];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [416198 36741 23952];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [419043 37293 24118];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [422140 38269 24301];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [425481 39480 24478];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [429193 40800 24656];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [432798 42112 24840];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [436319 43475 25015];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [439882 44818 25222];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [443086 45641 25394];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [446448 46689 25589];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [449860 47650 25779];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [453637 48924 25986];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [457219 50094 26169];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [461044 51296 26380];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [464596 52765 26567];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [468119 53698 26746];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [471772 54849 26957];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [475674 56189 27192];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [479825 57606 27419];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [483844 59077 27658];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [488236 61048 27888];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [492378 62901 28098];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [496253 64010 28293];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [500075 65142 28544];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [504281 66344 28816];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [508389 67479 29070];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [513219 69039 29349];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [517835 70176 29605];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [522387 71607 29870];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [526490 72446 30123];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [530380 72605 30375];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [534631 72559 30712];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [539670 73960 31034];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [545286 75231 31346];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [550757 76433 31650];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [555891 78221 31985];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [562705 79494 32320];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [568896 81226 32616];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [574856 82653 32953];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [581824 84914 33299];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [588648 87017 33714];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [596941 90230 34113];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [604952 93174 34478];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [612772 95978 34864];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [620491 98502 35298];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [628780 101795 35738];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [637712 106079 36160];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [646164 109185 36579];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [654936 112664 36985];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [663800 116439 37409];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [673250 120265 37832];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [682486 123866 38291];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [692949 128559 38749];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [703288 133392 39202];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [715068 139299 39664];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [726585 144898 40121];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [738322 151098 40582];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [749525 155744 41034];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [762068 161757 41493];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [775121 168443 41979];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [788473 175238 42461];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [801894 181970 42941];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [815117 187996 43417];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [828377 195456 43896];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [841308 200845 44327];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [854361 206114 44802];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [866821 211160 45255];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [880542 217089 45738];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [894385 222572 46207];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [908346 228382 46689];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [922397 235237 47095];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [935799 239482 47486];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [948749 242583 47874];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [962070 245673 48246];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [975951 249360 48628];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [989572 252528 48990];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1003494 254831 49348];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [1016835 259034 49695];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1028986 259262 50016];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1040547 259439 50310];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1051374 257825 50594];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1062397 257256 50917];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1072620 255769 51212];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1083023 253360 51496];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1092407 252827 51727];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1100818 248016 51949];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1108269 243803 52196];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1115770 240092 52447];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1123474 237528 52670];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1131077 233764 52883];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1138530 228922 53095];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1145651 226904 53273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1152072 222681 53448];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1158384 219705 53625];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1164535 216353 53816];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1170743 212742 54003];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1177004 208654 54156];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1183182 204189 54308];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1189203 201027 54440];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1194963 197841 54574];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1200465 193912 54693];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1206373 190808 54814];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1212481 188127 54946];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1218753 184944 55095];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1225142 181086 55223];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1231429 180522 55337];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1237474 177710 55438];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1243434 174876 55540];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1249507 173120 55650];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1255620 170844 55748];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1261903 166552 55830];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1268263 161777 55933];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1274514 159842 56108];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1280438 156872 56108];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1286406 155348 56171];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1292614 154616 56262];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1299022 154197 56360];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1305339 154494 56457];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1311810 153673 56538];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1318295 154663 56621];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1324395 154454 56717];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1330411 154471 56803];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1336217 153832 56886];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1342134 153883 56973];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1348316 153447 57057];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1354520 152821 57150];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1360832 151951 57225];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;







end