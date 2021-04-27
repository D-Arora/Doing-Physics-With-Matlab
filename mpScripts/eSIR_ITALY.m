%eSIR_ITALY.m

clear
close all
clc


global a b


% Model Parameters
  cn = 'ITALY';
  a = 0.18;
  b = 0.035;
  f = 1.8e5;
  
  k0 = 1.6e-5;
  D0 = 3.5e4;
  
  k1 = 0.06e-5;
  D1 = 7.5e4;
  ds1 = 2400;
   
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
  R(1) = 149/f;
  
  t = linspace(0,tMax,num);
  h = t(2) - t(1);
  
  
% SOLVE ODEs with RK4  ================================================
for n = 1 : num-1
   
   if n == 1250; S(n) = 0.18; end 
   if n > 1680 && n < 2230; S(n) = 0.40; end
   if n > 2229 && n < 2700; S(n) = 0.58; end
   
   if n > 2700; S(n) = 0.16; end
   
   
   if n > 3000 && n < 3200; S(n) = 0.3; end
   
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
  
  
  D(ds1:end) = D(ds1) + D1.*(1 - exp(-k1.*(R(ds1:end)-R(ds1))));
  
  C = R - D;
  Itot = I + R;
   
  Re = (a/b).*S;
  Idot = a.*S.*I - b.*I;
  
  
  
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
   text(-30,-7e3,'22 Jan Feb 2020')
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
   text(-30,-0.3,'15 Feb 2020')
   
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
%    set(gca,'xtick',0:50:tMax)
%    
%    
% figure(2)  % =========================================================
%    set(gcf,'units','normalized');
%    set(gcf,'position',[0.35 0.10 0.35 0.5]);
%    set(gcf,'color','w');
%    FS = 14;
%    
% subplot(2,2,1)   
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
% 
%    
% figure(9)
%    set(gcf,'units','normalized');
%    set(gcf,'position',[0.68 0.10 0.30 0.4]);
%    set(gcf,'color','w');
%    FS = 14;
%    
%    xlim([0 120])
%    ylim([0 250])
%    hh = 250; dh = -20;
%    
%    txt = 'ITALY:   START DATE   26 February 2020';
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
%    z = z + Ndays;
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

T = [470 455 12];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [655 593 17];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [889 822 21];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [1128 1049 29];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1701 1577 41];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [2036 1835 52];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2502 2263 79];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [3089 2706 107];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [3858 3296 148];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [4636 3912 197];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [5883 5016 233];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7375 6387 366];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [9172 7985 463];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10149 8514 631];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [12462 10590 827];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [15113 12838 1016];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [17660 14955 1266];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21157 17750 1441];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [24747 20608 1809];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [27980 23073 2158];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [31506 26062 2503];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [35713 28710 2978];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [41035 33190 3405];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [47021 37860 4032];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [53578 42687 4825];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [59138 46638 5476];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [63927 50418 6077];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [69176 54030 6820];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [74386 57521 7503];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [80589 62013 8215];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [86498 66414 9134];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [92472 70065 10023];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [97689 73880 10779];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [101739 75528 11591 ];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [105792 77635 12428];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [110574 80572 13155];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [115242 83049 13915];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [119827 85388 14681];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [124632 88274 15362];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [128948 91246 15887];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [132547 93187 16523];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [135586 94067 17127];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [139422 95262 17669];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [139422 96877 18279];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [147577 98273 18279];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [152271 100269 19468];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [156363 102253 19899];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [159516 103616 20465];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [162488 104291 21067];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [165155 105418 21645];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [168941 106607 22170];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [172434 106962 22745];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [175925 107771 23227];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [178972 108257 23650];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [181228 108237 24114];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [183957 107699 24648];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [187327 107699 25085];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [189973 106848 25549];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [192994 106527 25969];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [195351 105847 26384];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [197675 106103 26644];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [199414 105813 26977];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [201505 105205 27359];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [203591 104657 27682];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [205463 101551 27967];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [207428 100943 28236];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [209328 100704 28710];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [210717 100179 28884];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [211938 99980 29079];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [213013 98467 29315];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [214457 91528 29684];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [215858 89624 29958];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [217185 87961 30201];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [218268 84842 30395];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [219070 83324 30560];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [219814 82488 30739];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [221216 81266 30911];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [222104 78457 31106];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [223096 76440 31368];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [223885 72070 31610];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [224760 70187 31763];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [225435 68351 31908];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [225886 66553 32007];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [226699 65129 32169 ];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [227364 62752 32330];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [228006 60960 32486];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [228658 59322 32616];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [229327 57752 32735];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [229858 56594 32785];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [230158 55300 32877];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [230555 52942 32955];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [231139 50966 33072];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [231732 47986 33142];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [232248 46175 33229];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [232664 43691 33340];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [232997 42075 33415];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [233197 41367 33475];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [233515 39893 33530];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [233836 39297 33601];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [234013 38429 33689];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [234531 36976 33774];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [234801 35877 33846];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [234998 35262 33899];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [235278 34730 33964];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [235561 32872 34043];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [235763 31710 34114];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [236142 30637 34167];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [236305 28997 34223];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [236651 27485 34301];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [236989 26274 34345];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [237290 25909 34371];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [237500 24569 34405];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [237828 23925 34448];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [238159 23101 34514];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [238011 22353 34561];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [238275 21212 34610];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [238499 20972 34634];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [238720 20637 34657];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [238833 19573 34614];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [239410 18655 34644];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [239706 18303 34678];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [239961 17638 34708];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [240136 16836 34716];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [240310 16681 34738];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [240436 16496 34744];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [240578 15563 34767];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [240760 15255 34788];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [240961 15060 34818];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [241184 14884 34837];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [241419 14621 34854];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [241611 14642 34861];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [241819 14709 34869];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [241956 14242 34899];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [242149 13595 34914];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [242363 13459 34926];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [242639 13428 34938];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [242827 13303 34945];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [243061 13179 34954];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [243230 13157 34967];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [243344 12919 34984];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [243506 12493 34997];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [243736 12473 35017];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [243969 12456 35028];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [244216 12368 35042];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [244434 12440 35045];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [244624 12404 35058];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [244752 12248 35073];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [245032 12322 35082];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [245338 12404 35092];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [245590 12301 35097];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [245864 12442 35102];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [246118 12565 35107];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [246286 12581 35112];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [246488 12609 35123];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [246776 12616 35129];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [247158 12230 35132];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [247537 12422 35141];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [247832 12457 35146];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [248070 12456 35154];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [248229 12474 35166];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [248419 12482 35171];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [248803 12646 35181];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [249204 12694 35187];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [249756 12924 35190];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [250103 12953 35203];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [250566 13263 35205];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [250825 13368 35209];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [251237 13561 35215];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [251713 13791 35225];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [252235 14081 35231];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [252809 14249 35234];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [253525 14406 35392];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [253915 14733 35396];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [254235 14867 35400];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [254636 15089 35405];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [255278 15360 35412];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [256118 16014 35418];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [257065 16678 35427];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [258136 17503 35430];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [259345 18438 35437];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [260298 19195 35441];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [261174 19714 35446];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [262540 20753 35458];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [263949 21932 35463];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [265409 23035 35472];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [266853 23156 35473];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [268218 24205 35477];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [269214 26078 35483];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [270189 26754 35491];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [271515 27817 35497];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [272912 28915 35507];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [274644 30099 35518];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [276338 31194 35534];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [277634 32078 35541];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [280153 33789 35563];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [281583 34734 35577];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [283180 35708 35587];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [284796 36767 35597];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [286297 37503 35603];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [287753 38509 35610];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [288761 39187 35624];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [289990 39712 35633 ];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [291442 40532 35645];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [293025 41413 35658];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [294932 42457 35668];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [296569 43161 35692];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [298743 44098 35707];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [299506 45079 35724];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [300897 45489 35738];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [302537 46114 35758];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [304323 46780 35781];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [306235 47718 35801];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [308104 48593 35818];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [309870 49618 35835];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [311364 50323 35851];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [313011 50630 35875];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [314861 51263 35894];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [317409 52294 35918];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [319908 53997 35941];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [322751 55566 35968];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [325329 57429 35986];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [327586 58903 36002];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [330263 60134 36030];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [333940 62576 36061];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [338398 65952 36083];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [343770 70110 36111];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [349494 74829 36140];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [354950 79075 36166];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [359569 79792 36205];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [365467 87193 36246];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [372799 92445 36289];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [381602 99266 36372];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [391611 107312 36427];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [402536 116935 36474];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [414241 126237 36543];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [423578 134003 36616];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [434449 142739 36705];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [449648 155442 36832];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [465726 169302 36968];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [484869 186002 37059];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [504509 203182 37210];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [525782 222241 37338];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [542789 236684 37479];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [564778 255090 37700];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [589766 276457 37905];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [616595 299191 38122];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [647674 325786 38321];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [679430 351386 38618];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [709335 378129 38826];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [731588 396512 39059];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [759829 418142 39412];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [790377 443235 39764];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [824879 472348 40192];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [862681 499118 40638];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [902490 532536 41063];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [935104 558636 41394];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [960373 573334 41750];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [995463 590110 42330];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1028424 613358 42953];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1066401 635054 43589];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1107303 663926 44139];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1144552 688435 44683];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1178529 712490 45229];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1205881 717784 45733];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1238072 733810 46464];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1272352 743168 47217];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1308528 761671 47870];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1345767 777176 48569];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1380531 791746 49261];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1408868 805947 49823];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1455022 798386 51306];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1480874 791697 52028];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1538301 787893 53677];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1564532 789308 54363];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1585178 795771 54904];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1601554 788471 55576];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1620901 779945 56361];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [1641610 761230 57045];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [1664829 759982 58038];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1688939 757702 58852];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [1709991 754169 59514];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [1728878 755306 60078];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1742557 748819 60606];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1757394 737525 61240];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1770149 710515 61739];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1787147 696527 62626];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1805873 690323 63387];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1825775 684848 64036];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1843712 686031 64520];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1855737 675109 65011];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1870576 667303 65857];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1888144 645706 66537];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1906377 635343 67220];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1906377 635343 67220];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1938083 620166 68447];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1953185 622760 68799];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1964054 613582 69214];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1977370 605955 69842];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1991278  598816 70395];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2009317 593632 70900];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2028354 579886 71359];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2038759 580941 71620];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2047696 581760 71925];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2056277 575221 72370];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2067487 568728 73029];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2083689 564395 73604];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2107166 569896 74159];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2129376 574767 74621];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2141201 577062 74985];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2155446 576214 75332];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2166244 570458 75680];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2181619 596161 76329];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2201945 568712 76877];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2220361 571055 77291];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2237890 570389 77911];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2257866 572842 78394];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2275491 579932 78755];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2289021 575779 79203];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2303263 570040 79819];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2319036 564774 80326];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2336279 561380 80848];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2352423 558068 81326];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2368733 557717 81800];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2381277 553374 82177];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2390101 547058 82554];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2400598 535524 83157];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2414156 523553 83681];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2412505 516568 84202];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2441854 502053 84674];
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