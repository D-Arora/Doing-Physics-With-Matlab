%eSIR_KOREA.m

clear
close all
clc


global a b

% Country and start date
  cn = 'South Korea';   
  DS = '26 Feb 2020';            % 
  zSTART = datenum([2020 2 26]);

% Model Parameters
  a = 0.42;
  b = 0.0405;
  f = 1.65e4;
  
  k0 = 0.9e-4;
  D0 = 425;
  
  k1 = 0.5e-4;
  D1 = 520;
  ns1 = 2000;
  
  k2 = 0.7e-4;
  D2 = 750;
  ns2 = 2950;
  
   
  num = 5000;
  tMax = 500;
  
  
  
% Initialize matrices
  S  = zeros(num,1);       % Susceptible population
  I = zeros(num,1);        % Active infected population
  R = zeros(num,1);        % Removals form infected population  
 % C = zeros(num,1);        % ReCovered population
 % D = zeros(num,1);        % Dead population
  
  
  S(1) = 7;
  I(1) = 12e-2;
  R(1) = 36/f;
  t = linspace(0,tMax,num);
  h = t(2) - t(1);
  
  
% SOLVE ODEs with RK4  ================================================
for n = 1 : num-1
    
%    if n == 960; S(n) = 0.12; end  
%    if n == 1020; S(n) = 0.26; end
    if n > 700  &&  n < 1501; S(n) = 0.5; end 
    if n > 1600 &&  n < 1701; S(n) = 5; end 
    if n > 2200 &&  n < 2401; S(n) = 1.3; end 
    if n > 2600 &&  n < 2801; S(n) = 2.6; end
    if n > 2800 &&  n < 3101; S(n) = 2.2; end
    
    if n > 3100; S(n) = 0.55; end
        
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
  
  D(ns1:ns2) = D(ns1) + D1.*(1 - exp(-k1.*(R(ns1:ns2)-R(ns1))));
  
  D(ns2:end) = D(ns2) + D2.*(1 - exp(-k2.*(R(ns2:end)-R(ns2))));
  
  C = R - D;
  Itot = I + R;
  
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
   FS = 16;
 
  subplot(5,1,1)
   x1 = 1:Ndays; y1 = Idtot;
   x2 = t;       y2 = Itot;
   yL =  'Total  infections';
   plot1(x1,y1,x2,y2,tMax,cn,DS,FS, yL)
   title(cn)
  
   
  subplot(5,1,2)
   x1 = 1:Ndays; y1 = Id;
   x2 = t;       y2 = I;
   yL =  'Active  infections';
   plot1(x1,y1,x2,y2,tMax,cn,DS,FS, yL)
   
  subplot(5,1,3)
   x1 = 1:Ndays; y1 = Rd;
   x2 = t;       y2 = R;
   yL =  'Removed';
   plot1(x1,y1,x2,y2,tMax,cn,DS,FS, yL)
  
  subplot(5,1,4)
   x1 = 1:Ndays; y1 = Cd;
   x2 = t;       y2 = C;
   yL = 'reCoveries';
   plot1(x1,y1,x2,y2,tMax,cn,DS,FS, yL)
     
   
  subplot(5,1,5) 
   x1 = 1:Ndays; y1 = Dd;
   x2 = t;       y2 = D;
   yL = 'Deaths';
   plot1(x1,y1,x2,y2,tMax,cn,DS,FS, yL) 
   xlabel('Days elapsed')
   text(-30,-500,DS)
 
       
   
figure(2)  % =========================================================
   set(gcf,'units','normalized');
   set(gcf,'position',[0.52 0.01 0.40 0.85]);
   set(gcf,'color','w');
   FS = 16;
   
   zz = find(t >Ndays,1);

%   subplot(5,1,5) 
%    x1 = t; y1 = Dd;
%    x2 = t;       y2 = D;
%    yL = 'Deaths';
%    plot1(x1,y1,x2,y2,tMax,cn,DS,FS, yL) 
%    xlabel('Days elapsed')
%    text(-30,-500,DS)  
   
subplot(5,1,1)   
   plot(t,S,'r','linewidth',2)
   xlabel('Days elapsed')
   ylabel('R_E');
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

   Htext = plot(t(zz),S(zz),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
   text(-30,-3,'15 Feb 2020')
   
subplot(5,1,2)
   plot(t,Idot,'m')
   hold on
   plot(t,gradient(I,h),'r','linewidth',2)
   ylabel('dI/dt')
   xlabel('Days elapsed')
   
   Htext = plot(t(zz),b.*(S(zz)-1).*I(zz),'ro');
   set(Htext,'markersize',6,'MarkerFaceColor','r')
      
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   grid on
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')

 subplot(5,1,3)   
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
   
subplot(5,1,4)   
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
   
  subplot(5,1,5)   
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
   
   hh = hh+2*dh; 
   txt = 'DATA';
   text(0,hh,txt,'fontsize',12)
   
   z = zSTART;
   z = z + Ndays+1;
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
     txt = sprintf('percent: Active = %2.1f  ReCoveries = %2.1f    Deaths = %2.1f   ',pI,pC,pD);
     text(2,hh,txt,'fontsize',12)   
        
     axis off
    


% FUNCTIONS  ===========================================================

function graph1 = plot1(x1,y1,x2,y2,tMax,cn,DS,FS,yL)

   plot(x1,y1,'b+','linewidth',1);
   hold on
   plot(x2,y2,'r','linewidth',2)
   
   xlim([0 tMax])
   yMax = ylim;
   z = yMax(2);
   plotMonths(z)
   
   ylabel(yL)
   set(gca,'xtick',0:50:tMax)
   set(gca,'fontsize',FS)
   set(gca,'fontName','times')
   grid on; box on;
   
   

end

function FS = SDOT(t,S,I,R)
  global a
    FS = -a*S*I;
end


function FI = IDOT(t,S,I,R)
  global a b
    FI = b*(S-1)*I;
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



% FUNCTIONS ============================================================
  
function covid = sir
c = 1;
T = [1261 1225 12];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [1766 1729 13];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [2337 2297 16];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [3150 3109 17];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [3736 3685 21];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [4335 4277 28];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [5186 5120 32];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [5621 5498 35];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [6284 6107 42];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [6593 6415 43];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7041 6875 48];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7313 7097 50];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7478 7178 53];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7513 7165 60];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7755 7362 60];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7869 7293 66];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [7979 7198 67];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8086 7180 72];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8162 7253 75];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8236 7024 75];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8320 6838 81];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8413 6789 84];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8565 6257 91];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8652 6325 94];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8799 6085 102];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8897 5884 104];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [8961 5684 111];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9031 5410 120];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9241 5281 126];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9241 4966 131];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;   
T = [9332 4665 139];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9478 4523 144];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9583 4398 152];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9661 4275 158];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9786 4216 162];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9887 4155 165];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [9976 3979 169];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10062 3867 174];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10156 3654 177];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10237 3591 183];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10284 3500  186];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10331 3445 192];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10384 3408 200];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10423 3246 204];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10450 3125 208];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10480 3026 211];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10512 2930 214];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10537 2873 217];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10564 2808 222];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10591 2750 225];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10613 2627 229];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10635 2576 230];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10653 2484 232];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10661 2385 234];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10674 2324 236];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10683 2233 237];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10694 2179 238];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10702 2051 240];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10708 1967 240];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
   T = [10718 1843 240];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10728 1769 242];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10738 1731 243];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10752 1654 244];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10761 1593 246];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10765 1459 247];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10774 1454 248];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10780 1407 250];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10793 1360 250];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10801 1332 252 252];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10804 1267 254];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10806 1218 255];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10870 1135 256];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10822 1082 256];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10840 1016 256];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10874 1008 256];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10909 1021 256];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10936 1008 258];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10962 1008 259];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [10991 969 260];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11018 937 260];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11037 924 262];
  covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11065 898 263];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11065 898 263];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11037 924 262];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11122 723 264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11122 723 264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11142 716 264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11165 705 266];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11190 711 266];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11206 713 267];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11225 681 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11265 701 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11344 735 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11402 770 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11441 774 269];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11468 793 270];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11503 810 271];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11541 823 272];
   covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11590 850 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11629 857 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1; 
T = [11668 889 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11719 915 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11776 951 273];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11852 989 274];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11902 1015 276];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11902 1017 276];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [11947 1017 276];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12003 1057 277];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12051 1083 277];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12085 1090 277];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12121 1114 277];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12155 1117 278];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12294 1145 279];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12257 1177 280];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12306 1191 280];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12373 1237 280];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12421 1273 280];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12484 1295 281];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12535 1324 281];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12563 1307 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12602 1148 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12653 1054 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12715 1069 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12757 1046 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12800 981 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12850 955 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12904 938 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [12967 926 282];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13030 936 283];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13091 976 283];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13137 1005 284];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13181 982 285];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13244 989 285];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13293 987 287];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13338 985 288];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13373 941 288];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13417 950 289];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13479 986 289];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13512 941 289];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13551 914 289];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13612 925 291];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13672 919 293];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13711 898 294];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13745 894 295];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13771 903 296];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;% T = [];
T = [13816 877 296];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13879 884 297];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13938 883 297];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [13979  864 298];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14092 928 298];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14150 962 298];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14175 971 299];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14203 896 300];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14251 882 300];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14269 837 300];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14305 821 301];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14336 802 301];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14366 806 301];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14389 808 301];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14423 770 301];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14456 748 302];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14499 696 302];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14519 673 303];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14562 629 304];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14598 651 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14626 663 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14660 626 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14714 623 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14770 648 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [14873 705 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [15039 833 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [15318 1103 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [15515 1293 305];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [15761 1521 306];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [16058 1746 306];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [16346 1976 307];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [16670 2241 309];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [17002 2524 309];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [17399 2890 309];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [17665 3137 309];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [17945 3349 310];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [18265 3585 312];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [18706 3932 313];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [19077 4210 316];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [19400 4314 323];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [19699 4473 323];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [19947 4650 324];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [20182 4660 324];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [20449 4767 326];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [20644 4786 329];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [20842 4728 331];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21010 4668 333];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21177 4697 334];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21296 4663 336];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21432 4455 341];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21588 4221 344];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21743 4037 346];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [21919 3953 350];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22055 3671 355];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22176 3592 358];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22285 3433 363];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22391 3146 367];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22504 2827 367];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22657 2742 372];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22783 2625 377];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22893 2545 378];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [22975 2434 383];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23045 2412 385];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23106 2277 388];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23216 2178 388];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23341 2116 393];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23455 2082 395];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23516 1951 399];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23611 1962 401];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23661 1963 406];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23699 1822 407];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23812 1809 413];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23889 1808 415];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [23952 1803 416];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24027 1820 420];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24091 1825 421];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24164 1856 422];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24239 1734 422];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24353 1594 425];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24422 1532 427];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24476 1479 428];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24548 1494 430];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24606 1481 432];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24703 1541 433];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24805 1508 434];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24889 1421 438];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [24988 1467 439];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25035 1414 441];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25108 1407 443];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25199 1443 444];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25275 1463 444];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25333 1420 447];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25424 1390 450];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25543 1443 453];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25698 1526 455];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25775 1484 457];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25836 1510 457];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [25955 1593 457];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26043 1602 460];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26146 1621 461];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26271 1641 462];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26385 1695 463];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26511 1736 464];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26732 1869 468];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26732 1869 468];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26807 1825 472];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [26925 1835 474];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27050 1840 475];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27195 1898 476];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27284 1981 477];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27427 1981 478];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27553 2044 480];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27653 2008 485];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27799 2046 487];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [27942 2051 487];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [28133 2108 488];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [28338 2210 492];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [28546 2362 493];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [28769 2516 494];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [28392 2644 494];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [29311 2842 496];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [29654 3058 498];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [30017 3253 501];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [30403 3535 503];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [30733 3762 505];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [31004 3956 509];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [31353 4121 510];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [31735 4397 513];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [32318 4853 515];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [32887 5268 516];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [33375 5504 522];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [33824 5769 523];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [34201 6022 526];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [34652 6241 526];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [35163 6572 526];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [36332 7185 536];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [36332 7185 536];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
 T = [36915 7458 540];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [38161 8311 549];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [38755 8558 552];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [39432 8699 556];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [40098 8897 564];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [40786 9057 572];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [41736 9665 578];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [42766 9638 580];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [43484 10795 587];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [44364 11205 600];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [45442 11883 612];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [47515 12888 634];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [47515 12888 645];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [48750 13577 659];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [49665 14269 674];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [50591 14738 698];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [51460 14810 722];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [52550 15085 739];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [53533 15352 756];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [54770 15949 773];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [56872 17024 808];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [56872 17024 808];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [57680 17593 819];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [58725 17163 859];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [59773 17459 879];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [60740 17569 900];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [61769 17899 917];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [62593 18073 942];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [63244 17775 962];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [64264 18043 981];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [64979 17800 1007];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [65818 17796 1027];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [66686 17991 1046];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [67358 17908 1081];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [67999 17575 1100];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [68664 17130 1125];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [69114 15422 1140];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [69651 14917 1165];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [70212 14391 1185];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [70728 13761 1195];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [71241 13488 1217];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [71820 13090 1236];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [72340 12838 1249];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [72729 12742 1264];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [73115 12364 1283];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [73518 12038 1300];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [73918 11756 1316];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
T = [74262 11519 1328];
    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;
% T = [];
%    covid(c,1) = T(1); covid(c,2) = T(2); covid(c,3) = T(3); c = c+1;










end