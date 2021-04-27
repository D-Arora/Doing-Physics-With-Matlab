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
  a = 0.08;
  b = 0.032;
  f = 35e5;
  
  k0 = 0.07e-5;
  D0 = 75e3;
  
  tMax = 300;
  
  I(1) = 50e-5;
  R(1) = 100e-4;
  
% Starting Date
  zSTART = datenum([2020 3 14]);  
  
  
% Setup  ==============================================================
  S(1) = 1;
  t = linspace(0,tMax,num);
  h = t(2) - t(1);
  
  
% SOLVE ODEs with RK4  ================================================
for n = 1 : num-1
    
%   if n == 2200; S(n) = 0.7; end
%   if n == 2300; S(n) = 0.7; end
%   if n == 2400; S(n) = 0.7; end
%   
%  if n == 1733; S(n) = 0.7; end
%  if n == 1800; S(n) = 0.7; end
%  if n == 1867; S(n) = 0.8; end
%  if n == 1933; S(n) = 0.8; end
%  if n == 2000; S(n) = 0.8; end  
%  if n == 2067; S(n) = 0.8; end  
%  if n == 2133; S(n) = 0.8; end  
%  if n == 2400; S(n) = 0.7; end  
%  if n == 2200; S(n) = 0.7; end  
%  if n == 2000; S(n) = 0.7; end  
 
%   if n > 2200 && n < 2500
%            S(n) = 0.60;
%   end 

%     if n > 2199 && n < 2500
%          S(n) = 0.70;
%     end 
    
%  if n > 2300 && n < 2600 %2400  %2834
%      S(n) = 0.30;
%  end

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
  
% D = 0.055.*R;
  
  C = R - D;
  Itot = I + R;
  
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
  
  
% GRAPHICS  ===========================================================

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.0 0.10 0.35 0.5]);
   set(gcf,'color','w');
   FS = 14;
   
subplot(2,2,1)   
   plot(1:Ndays,Idtot,'b+','linewidth',1)
   xlabel('Days elapsed')
   ylabel('Total  infections')
   title(cn)
   hold on
   plot(t,Itot,'r','linewidth',1)
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
 
 subplot(2,2,2) 
   plot(1:Ndays,Id,'b+','linewidth',1)
   grid on; box on;
   xlabel('Days elapsed')
   ylabel('Active  infections')
   hold on
   plot(t,I,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
  
  subplot(2,2,3)   
   plot(1:Ndays,Cd,'b+','linewidth',1)
   grid on; box on;
   xlabel('Days elapsed')
   ylabel('Recoveries')
   hold on
   plot(t,C,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   
  subplot(2,2,4)   
   plot(1:Ndays,Dd/1000,'b+','linewidth',1)
   grid on; box on;
   xlabel('Days elapsed')
   ylabel('Deaths/1000')
   hold on
   plot(t,D/1000,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   
   
figure(2)  % =========================================================
   set(gcf,'units','normalized');
   set(gcf,'position',[0.35 0.10 0.35 0.5]);
   set(gcf,'color','w');
   FS = 14;


   subplot(2,2,1)   
   plot(t(1:4880),Re(1:4880),'r','linewidth',2)
   xlabel('Days elapsed')
 %  ylabel('Susceptible  S')
   ylabel('Reproduction  No.   R_e')
   title(cn)
   xlim([0 tMax])
 %  ylim([0 1.1])
   ylim([0 1.1*max(Re)])
 %  set(gca,'ytick',0:0.2:1)
   set(gca,'xtick',0:50:tMax)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   hold on
   grid on; box on
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   
 %  plot([0 tMax],[b/a b/a],'m','linewidth',1)
   plot([0 tMax],[1 1],'m','linewidth',1)
  % txt = sprintf('S_C = %2.2f',b/a);
  % Htext = text(110,1.6*b/a,txt);
  % set(Htext,'fontName','times','fontSize',12,'BackgroundColor','w')
   
   z = find(t >Ndays,1);
   Htext = plot(t(z),Re(z),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
   
% subplot(2,2,1)   
%    plot(t,S,'r','linewidth',2)
%    xlabel('Days elapsed')
%    ylabel('Susceptible  S')
%    title(cn)
%    xlim([0 tMax])
%    ylim([0 1.1])
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
%    plot([0 tMax],[b/a b/a],'m','linewidth',1)
%    txt = sprintf('S_C = %2.2f',b/a);
%    Htext = text(30,0.5*b/a,txt);
%    set(Htext,'fontName','times','fontSize',12,'BackgroundColor','w')
%    
%    z = find(t >Ndays,1);
%    Htext = plot(t(z),S(z),'ro');
%    set(Htext,'markersize',6,'MarkerFaceColor','r')
%    
 subplot(2,2,2)   
   plot(Rd,Dd,'b+','linewidth',1)
   grid on; box on;
   xlabel('Removals')
   ylabel('Deaths')
   hold on
   plot(R,D,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   
subplot(2,2,3)   
   plot(Id,Rd,'b+','linewidth',1)
   grid on; box on;
   xlabel('Active infections')
   ylabel('Removals')
   hold on
   plot(I,R,'r','linewidth',2)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   
  subplot(2,2,4)   
   plot(Idtot,Rd,'b+','linewidth',1)
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
   Htext = text(20,hh,txt,'fontsize',12);
   set(Htext,'color','k')
   z = zSTART;
   z = datetime(z,'ConvertFrom','datenum');
   z.Format = 'dd-MMM-yyyy';
   txt = cellstr(z) ; 
   text(45,hh,txt,'fontsize',12)
   
   hh = hh+2*dh;
   txt = 'MODEL PARAMETERS';
   text(0,hh,txt,'fontsize',12)
   
    hh = hh+dh;
   txt = sprintf('I(0) = %3.2e  R(0) = %2.2e  f = %2.2e  a = %2.3f  b = %2.3f', I(1)/f, R(1)/f, f,a ,b);
   text(0,hh,txt,'fontsize',12)
   
   hh = hh+dh;
   txt = 'D = D_0 (1-exp(-k_0 R))';
   text(0,hh,txt,'fontsize',12)
   txt = sprintf('D_0 = %3.2e   k_0 = %3.2e', D0, k0);
   text(50,hh,txt,'fontsize',12)
   
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
       m2 = [2020 04  1];    % Apr
       m3 = [2020 05  1];    % May
       m4 = [2020 06  1];    % Jun
       m5 = [2020 07  1];    % Jul
       m6 = [2020 08  1];    % Aug
       m7 = [2020 09  1];    % Sep
       m8 = [2020 10  1];    % Oct
       m9 = [2020 11  1];
       m10 = [2020 12  1];
       m11 = [2021 1  1];
       
   
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
      
       for n = 2:11
        plot([tm(n), tm(n)], [0, yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
        hold on
      end
      
      ys = 0.94;
      M = 'AMJJASOND';
      for c = 1:length(M)
        X = 26+(c-1)*30;  
        text(X, ys*yMax,M(c))
      end

      
        
%       for n = 2:11
%         plot([tm(n), tm(n)], [0, yMax],'linewidth',1,'color', [0.5 0.5 0.5]);
%         hold on
%       end
%       
%       ys = 0.9;
%       text(24, ys*yMax,'A')
%       text(54, ys*yMax,'M')
%       text(84, ys*yMax,'J')
%       text(114,ys*yMax,'J')
%       text(144,ys*yMax,'A')
%       text(174,ys*yMax,'S')
      
   end

