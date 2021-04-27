% eSIR.m

clear
close all
clc


global a b

% Initialize matrices
  num = 5000;
  S  = zeros(num,1);       % Effective Reproduction number
  I = zeros(num,1);        % Active infected population
  R = zeros(num,1);        % Removals form infected population 


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

flagC = 11;

% COUNTRIES  ==========================================================   
  switch  flagC

      
%   TEXAS *************************************************************
    case 3
         cn = "TEXAS:    1 Apr 20  to  13 Aug 21";
         a = 0.130;
         b = 0.040;
         f = 72e4;
         
         S(1) = 1.00;
         I(1) = 8*1.8e-5;
         R(1) = 5*1.6e-6;
  
         k0 = 1.10e-5; D0 = 5e3;
         k1 = 0.10e-5; D1 = 28000; ds1 = 1300;
         k2 = 0.20e-5; D2 = 15000; ds2 = 2500;
          
         Ndays  = find(covid == 0,1) - 1;  
         Idtot  = covid(1:Ndays,3*flagC -2);
         Id     = covid(1:Ndays,3*flagC - 1);
         Dd     = covid(1:Ndays,3*flagC); 
         
         Cd = Idtot - Id - Dd;
         Rd = Cd + Dd;
         
              
  end
  
 
 % SOLVE ODE  ========================================================= 
for n = 1 : num-1
        
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
   
   S(n+1) = ERN(flagC,S(n+1),n) ; 
end
  
  I = f.*I;
  R = f.*R;

  D = D0.*(1 - exp(-k0.*R));
  D(ds1:ds2) = D(ds1) + D1.*(1 - exp(-k1.*(R(ds1:ds2)-R(ds1))));
  D(ds2:end) = D(ds2) + D2.*(1 - exp(-k2.*(R(ds2:end)-R(ds2))));
  

  C = R - D;
  Itot = I + R;
  
  Idot = a.*S.*I - b.*I;
  
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
  % plot(t,C,'r','linewidth',2)
   yLabel = 'R reCoveries';
   graph(tMax, FS, yLabel) 
   
   %%
   mx = 8600;
   xx = 100:150;
   bx = Rd(133) - mx*133;
   yy = mx.*xx + bx;
   
   plot(xx,yy)
   
   
   %%
   
  
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
   plot(t,S,'r','linewidth',2)
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

   Htext = plot(t(zz),S(zz),'ro');
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



function z = ERN(flagC,V,n)
      z = V;
switch flagC
    
% TEXAS *************************************************************** 
    case 3 
%        if n > 500  && n < 1001; Re = 2.0; end
%        if n > 1001 && n < 1201; Re = 1.5; end
%        if n > 1200 && n < 1901; Re = 0.8;   end
%        if n > 1900 && n < 2201; Re = 1.3; end   
%        if n > 2200 && n < 2401; Re = 1.6; end
%        if n > 2400 && n < 2601; Re = 1.6; end
%        if n > 2600 && n < 2701; Re = 1.6; end
%        if n > 2700 && n < 3001; Re = 1.4; end
% %      if n > 3000 && n < 3021; Re = 0.60*2.6; end

     
    if n > 1000 && n < 1101; z = 0.60; end
    if n > 1300 && n < 1501; z = 0.280; end
%     if n > 2380 && n < 2651; z = 0.56; end
%     if n > 2650 && n < 3001; z = 0.46; end
%     if n > 3000 && n < 3021; z = 0.60; end
%     

end     
     
end




