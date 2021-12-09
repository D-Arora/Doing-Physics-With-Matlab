% covid21ReB.m


% ******  https://ourworldindata.org/covid-vaccinations


% AUSTRALIA DATA
%   **** https://covidlive.com.au/report/daily-vaccinations/aus
%        https://covidlive.com.au/cases
%    New Daily Cases Smoothed  G / Total Vaccinations  AI
%   covid21H   DC: UK col 3 /  AUS  col6   
%    https://ourworldindata.org/coronavirus/country/australia

% Australia Hospitalizations      covid21H     col 3 
%     

% Australia Total Vaccinations    covid21H   col 4
%
   

% UK  Hospitalizations     covid21H  col 1
%     https://ourworldindata.org/grapher/uk-daily-covid-admissions
%    https://coronavirus.data.gov.uk/details/healthcare

%  https://coronavirus.data.gov.uk/

% UK  Vaccinations        covid21H   col 2
%     https://ourworldindata.org/grapher/cumulative-covid-vaccinations?tab=chart&country=~GBR


% https://ourworldindata.org/coronavirus/country/australia
% C country / D date / G new cases smoothed / H total deaths / R icu patients
% T hosiptal patients / AA total test / AI total-vaccinations / AJ people
% vacinated / AU population


% covidMinModel.docx
% covid210909.docx

% covid19.xlsx


% ========================================================================
clear 
close all
clc


tic

% Select COUNTRY  ====================================================
% AUSTRALIA 1  /  UK 2

  flagC = 2;


% SETUP  ==============================================================
  nT = 5000; tMax = 500;
  t = linspace(0,tMax,nT);
  h = t(2)-t(1);

  I = zeros(nT,1);
  R = zeros(nT,1);
  RE = zeros(nT,1);
  delta = zeros(nT,1);
  b = zeros(nT,1);


% LOAD DATA from covid21.mat covid21 ================================== 
   % load('covid21.mat');
   % load('covid21H.mat');
     load('covid21M.mat')
    
   switch flagC 
    case 1
      cn = '      AUSTRALIA';
      pop = 25862611;  
      Ndays  = find(covid21M(:,1) == 0,1) - 1;  
      Idtot  = covid21M(1:Ndays,1);
      Id     = covid21M(1:Ndays,2);
      Dd     = covid21M(1:Ndays,3); 
      Hd     = covid21M(:,4);    % Data: hosipitalizations
      Vd     = covid21M(:,5);    % Data: Number of vaccinations 
      DCd    = covid21M(:,6);   % Data: Daily new cases of infection
       
     case 2  
      cn = '    U.K.';   
      pop = 68323609;
      Ndays  = find(covid21M(:,7) == 0,1) - 1;  
      Idtot  = covid21M(1:Ndays,7);
      Id     = covid21M(1:Ndays,8);
      Dd     = covid21M(1:Ndays,9); 
      Hd     = covid21M(:,10);    % Data: hosipitalizations
      Vd     = covid21M(:,11);    % Data: Number of vaccinations 
      DCd    = covid21M(:,12);   % Data: Daily new cases of infection
   end
    Cd = Idtot - Id - Dd;
    Rd = Cd + Dd;
    dayR = 1: Ndays;   
 
    I(1) = Id(1);
    R(1) = Rd(1);
  
% Estimate for b
    sRd = smooth(dayR,Rd,5);      %  Smooth Rd
    sId = smooth(dayR,Id,5); 
   for c = 1 : floor(Ndays/30)
      b((c-1)*300+1:c*300) = (sRd(c*30)-sRd((c-1)*30+1))/(30*Id((c-1)*30+15));
   end
      b(c*300:end) = b(c*300);
    
  bs = smooth(t,b,200);
  
  b = bs;
  
  if flagC == 4; b = 0.07.*ones(nT,1); end
%   %  bs = zeros(nT,1);
%     bs(1:Ndays) = gradient(sRd,1)./Id;
%     bs(Ndays:end) = bs(Ndays);
%   bs = smooth(dayR,bs,20);
%     
%     figure(77)
%     plot(t,bs)
%     
    % b0 = mean(b);
%       b0 = 0.02;
%       b = b0.*ones(nT,1);
      
% COUNTRY DEATHS
     D = Dd(1).*ones(nT,1);
     q = 1;
     
switch flagC  
     
 % 44444444444444444444444444444444444444444444444444444444444444444444444
  case 1    % AUSTRALIA
    RE(1) = 0.8;   
    a = 1e-5.*ones(nT,1);    
    b = 0.07.*ones(nT,1);
    delta(1000:1200) = 0.3e-2; 
    delta(1200:1400) = 0.2e-2;
    delta(1400:1600) = 0.1e-2;
    delta(1600:1800) = 0.1e-2;
    delta(1800:2000) = 0.3e-2;
    delta(2000:2200) = 0.3e-2;
    delta(2200:2600) = 1.5e-2;
    delta(2600:2800) = 3.0e-2;
    delta(2800:3000) = 3.3e-2;
    delta(3000:3200) = 1.5e-2;
    delta(3200:3700) = 5.0e-3;
    delta(3700:4300) = 30.0e-3;
    delta(4300:5000) = 10e-3;



     
    m(q) = 1; m(q+1) = 1100; d(q) = 10; k(q) = 1e-4;
    q = q+1;  m(q+1) = 1920; d(q) = 5; k(q) = 1e-4;
    q = q+1;  m(q+1) = 2000; d(q) = 20; k(q) = 1e-3;
    q = q+1;  m(q+1) = 2400; d(q) = 400; k(q) = 0.3e-4;
    q = q+1;  m(q+1) = 2500; d(q) = 400; k(q) = 0.2e-4;
    q = q+1;  m(q+1) = 2800; d(q) = 500; k(q) = 0.2e-4;
    q = q+1;  m(q+1) = 3200; d(q) = 500; k(q) = 0.3e-4;
    
    q = q+1;  m(q+1) = nT; d(q) = 20; k(q) = 1e-4;
    segments = length(k);
  
% 55555555555555555555555555555555555555555555555555555555555555555  
  case 2   % UK
    RE(1) = 6;
    a = 1e-7.*ones(nT,1);      
    
   %  b0 = 0.04;
   %   b = b0.*ones(nT,1);
     delta(50:200) = 5.5e-2;
     delta(200:1000) = 0.1e-2;
     delta(1400:1600) = 1.5e-2;
     delta(1600:1800) = 1.5e-2;
     delta(1800:2000) = 0.8e-2;
     delta(2000:2200) = 0.1e-2;
     delta(2200:2400) = 1e-2;
     delta(2400:2600) = 0.85e-2;
     delta(2600:2700) = 1.2e-2;
     delta(2700:2800) = 1e-2;
     delta(2800:2900) = 1.5e-2;
     delta(2900:3000) = 0.5e-2;
     delta(3000:3200) = 1.10e-2;
     delta(3200:3400) = 1.20e-2;

%     delta(3500:3800) = 3.2e-2;
%     delta(3800:4200) = 1.5e-2;
%     delta(4200:4400) = 0.8e-2;
%     delta(4400:5000) = 0.5e-2;
    
    D(1) = Dd(1);
    m(q) = 1; m(q+1) = 1000; d(q) = 52e3; k(q) = 2.2e-6;
    q = q+1;  m(q+1) = 2000; d(q) = 8e3; k(q) = 1e-6;
    q = q+1;  m(q+1) = 2500; d(q) = 8e3; k(q) = 1e-6;
    q = q+1;  m(q+1) = 3000; d(q) = 5e3; k(q) = 0.8e-6;
    q = q+1;  m(q+1) = 3500; d(q) = 20e3; k(q) = 0.5e-6;
    q = q+1;  m(q+1) = 3800; d(q) = 10e3; k(q) = 0.5e-6;

    q = q+1;  m(q+1) = nT; d(q) = 1e3; k(q) = 1e-6;
    segments = length(k);
end

 
% *******************************************************************
% CALCULATIONS  
% *******************************************************************
   delta = smooth(delta,200);
for c = 1 : nT-1
    
%    if c == 2530 && flagC == 4
%         I(c) = 16743;
%         R(c) = 52093 + 1076;
%    end     
    
   RE(c+1) = RE(c) - h*(a(c)*RE(c)*I(c)) + delta(c);
   I(c+1) = I(c) + h*(b(c)*(RE(c+1)-1)*I(c));
   R(c+1) = R(c) + h*b(c)*I(c);
end

% Deaths
  for c = 1 : segments
      D(m(c): m(c+1)) = D(m(c)) + d(c).*(1 - exp(-k(c).*(R(m(c) : m(c+1)) - R(m(c)))));
  end  
      Df = mean(Dd./Idtot);   % Death factor
  %    Df = 0.01;
% Total Infections
  Itot = I + R;
% ReCoveries      
  C = R - D;


% Deaths  
  Dr = Df.*Itot;                     % Predicted deaths
  Ddp = 100*Dd(Ndays)/Idtot(Ndays);  % Data: % deaths
  Dp = 100*D(end)/Itot(end);         % Model % deaths 
% Populations
  IP = 100*Itot(end)/pop;
% dates
  tS = datetime(2021,1,1);
  tC = tS + Ndays-1;

  

%  ******************************************************************* 
% HOSPITALIZATIONS / VACCINATIONS / DAILY NEW CASES 
%  ******************************************************************* 
   
%   HOSPITALIZATIONS
     Hs = Hd./max(Hd);    % Data: scaled hospitalizations
     maxI = max([max(Id), max(I)]);
     Ids = Id./maxI;   % Data: scaled current infections
     Is =  I./maxI;    % Model: scaled curremt infections
     ts1 = find(Hd == 0, 1) - 2; % Data:  time span
     if ts1 > Ndays; ts1 = Ndays; end
     ts = 1: ts1;         
    
%  Model prediction: Hospitalizations;
     tp1 = 10*(Ndays-2);         % Model: time span
     tp = tp1:nT; 
     Isp = I(tp)./maxI;
     sf = mean(Hd(ts1-7:ts1)./Id(ts1-7:ts1));
     Hp = I(tp).* sf;
     Hps =  Hp./max(Hd);

% VACCINATIONS    
    tv1 = find(Vd == 0, 1) - 1;
    tv = 1:tv1;
    Vs = Vd/(2*pop);
    
% Daily Cases
   Its = smooth(Itot,7);
   tDC = 1:500;
   DC = zeros(500,1);
   for c = 2:500
     DC(c) = Its(10*c) - Its(10*((c-1)));  
   end
   DC(1) = DC(2);
   DC = smooth(DC,7);
  

% DEATH RATE  ====================================================
  Ds = smooth(t,D,7);
  dDdt = 10*gradient(t,D,h);
  dDdt = smooth(t,dDdt,7);

  Davg = mean(Dd(Ndays-7:Ndays));
  Iavg = mean(Idtot(Ndays-7:Ndays));
  Dindex = 1*Davg/Iavg;
  D7 = Dindex.*Itot;


%  ******************************************************************* 
%  GRAPHICS 
%  *******************************************************************

  Z1 = [31	59	90	121	152	182	213	244	275	306	336	367	398	426	457	488];
  Z2 = 'JFMAMJJASONDJFMA';

  
% figure(1)  % 111111111111111111111111111111111111111111111111111111
%    set(gcf,'units','normalized');
%    set(gcf,'Position', [0.05 0.05 0.35 0.8500])
%    set(gcf,'color','w');
%    FS = 12;
%    width = 0.35; height = 0.20;
%    v3 = 0.15; v2 = 0.45; v1 = 0.75;
%    h1 = 0.10; h2 = 0.56;
%    
%    
% subplot('Position',[h1 v1 width height])
%   xtxt = 'days';
%   ytxt = 'I_{tot}';
%   xP = dayR; yP = Idtot; 
%   plot(xP,yP,'b+')
%   hold on
%   xP = t; yP = Itot; 
%   plot(xP,yP,'r','linewidth',2)
%  
%   xlim([0 500])
%   yMax = ylim;
%   ylim([0.6*yMax(1) yMax(2)])
%   
%   ax = gca;ax.FontSize = FS;
%   
%   Htitle = title(strcat(cn,{'     '},datestr(tC)));
%   set(Htitle,'fontweight','normal','fontsize',10)
%   plotMonths(xtxt, ytxt)
%   
%   
%   
% subplot('Position',[h2 v1 width height])
%   xtxt = 'days';
%   ytxt = 'I';
%   xP = dayR; yP = Id; 
%   plot(xP,yP,'b+')
%   hold on
%   xP = t; yP = I; 
%   plot(xP,yP,'r','linewidth',2)
%  
%   txt = sprintf('   I_{tot} / pop = %3.1f %% \n',IP);
%   xlim([0 500])
%   yMax = ylim;
%   ylim([0 yMax(2)])
%   
%   ax = gca;ax.FontSize = FS;
% 
%   plotMonths(xtxt, ytxt)
%   
%   title(txt,'fontweight','norma','fontsize',10)
% 
% 
% subplot('Position',[h1 v2 width height])
%   xtxt = 'days';
%   ytxt = 'R_C';
%   xP = dayR; yP = Cd; 
%   plot(xP,yP,'b+')
%   hold on
%   xP = t; yP = C; 
%   plot(xP,yP,'r','linewidth',2)
%  
%   xlim([0 500])
%   yMax = ylim;
%   ylim([0.6*yMax(1) yMax(2)])
%   
%   ax = gca;ax.FontSize = FS;
% 
%   plotMonths(xtxt, ytxt)
% 
%   
% subplot('Position',[h2 v2 width height]) 
%   xtxt = 'days';
%   ytxt = 'D / 1000';
%   xP = dayR; yP = Dd./1000; 
%   plot(xP,yP,'b+')
%   hold on
%   xP = t; yP = D./1000; 
%   plot(xP,yP,'r','linewidth',2)
% %      yP = Dr./1000;
% %      plot(xP,yP,'m','linewidth',2)
%    
%   xlim([0 500])
%   yMax = ylim;
%   ylim([0.6*yMax(2) yMax(2)])
%   
%   txt = sprintf('  D_{data} = %1.1f %%   D_{model} = %1.1f %% \n',Ddp, Dp);
% %  title(txt,'fontweight','normal','fontsize',10)
%   text(10,1.01*yMax(2),txt)
%   ax = gca;ax.FontSize = FS;
% 
%   plotMonths(xtxt, ytxt)
%   
% subplot('Position',[h1 v3 width height]) 
%   xtxt = 'I';
%   ytxt = 'R_M';grid on; box on;
%   xP = Id; yP = Rd; 
%   plot(xP,yP,'b+')
%   hold on
%   xP = I; yP = R; 
%   plot(xP,yP,'r','linewidth',2)
%  
%   yMax = ylim;
%   ylim([yMax(1) yMax(2)])
%   grid on; box on;
%   ax = gca;ax.FontSize = FS;
% 
%   xlabel(xtxt,'fontname','times','Fontsize',FS);
%   ylabel(ytxt,'fontname','times','Fontsize',FS);
%   
%   
% subplot('Position',[h2 v3 width height]) 
%   xtxt = 'R_M';
%   ytxt = 'D / 1000';grid on; box on;
%   xP = Rd; yP = Dd./1000; 
%   plot(xP,yP,'b+')
%   hold on
%   xP = R; yP = D./1000; 
%   plot(xP,yP,'r','linewidth',2)
%  
%   yMax = ylim;
%   ylim([yMax(1) yMax(2)])
%   grid on; box on;
%   ax = gca;ax.FontSize = FS;
% 
%   xlabel(xtxt,'fontname','times','Fontsize',FS);
%   ylabel(ytxt,'fontname','times','Fontsize',FS);
% 
%     
% figure(2)  % 222222222222222222222222222222222222222222222222222222
%    set(gcf,'units','normalized');
%    set(gcf,'Position', [0.5 0.05 0.35 0.8500])
%    set(gcf,'color','w');
%    FS = 12;
%    width = 0.35; height = 0.20;
%    v3 = 0.15; v2 = 0.45; v1 = 0.75;
%    h1 = 0.10; h2 = 0.56;
%    
%   fig = figure(2); 
%   left_color  = [0 0 1];
%   right_color = [0 0 0];
%   set(fig,'defaultAxesColorOrder',[left_color; right_color]);   
%    
% subplot('Position',[h1 v1 width height])
%   xtxt = 'days';
%   ytxt = 'R_E';
%   xP = t; yP = RE; 
%   plot(xP,yP,'r','linewidth',2)
%   hold on
%   plot([0 500],[1 1],'k','linewidth',1)
%   
%   xlim([0 500])
%   yMax = ylim;
%   ylim([0 yMax(2)])
%   
%   ax = gca;ax.FontSize = FS;
%   
%   Htitle = title(strcat(cn,{'     '},datestr(tC)));
%   set(Htitle,'fontweight','normal','fontsize',10)
%   plotMonths(xtxt, ytxt)
%    
%   
% subplot('Position',[h2 v1 width height])
%   xtxt = 'days';
%   ytxt = '\Delta';
%   xP = t; yP = delta; 
%   plot(xP,yP,'r','linewidth',2)
%  
%   xlim([0 500])
%   yMax = ylim;
%   ylim([0 yMax(2)])
%   
%   ax = gca;ax.FontSize = FS;
%   plotMonths(xtxt, ytxt)
%   
%    
% subplot('Position',[h1 v2 width height])
%   xtxt = 'days';
%   ytxt = 'a';
%   xP = t; yP = a; 
%   plot(xP,yP,'r','linewidth',2)
%  
%   xlim([0 500])
%   yMax = 1.1*max(a);
%   ylim([0 yMax])
%   
%   ax = gca;ax.FontSize = FS;
%   plotMonths(xtxt, ytxt)
% 
%   
% subplot('Position',[h2 v2 width height]) 
%   xtxt = 'days';
%   ytxt = 'b';
%   xP = t; yP = b; 
%   plot(xP,yP,'r','linewidth',2)
%  
%   xlim([0 500])
%   yMax = 1.1*max(b);
%   ylim([0 yMax])
%   
%   ax = gca;ax.FontSize = FS;
%   plotMonths(xtxt, ytxt)
% 
%   
% subplot('Position',[h1 v3 width height]) 
%     xtxt = 'days';
%     ytxt = 'scaled pop.';grid on; box on;
%    
%    hold on
%    xP = ts; yP = Ids(ts);
%       plot(xP,yP,'g-','linewidth',2)
%    xP = tp./10; yP = Isp;   
%       plot(xP,yP,'r-','linewidth',2)
%    xP = tp./10; yP = Hps;   
%       plot(xP,yP,'m-','linewidth',2)
%    xP = ts; yP = Hs(ts);
%       plot(xP, yP,'b-','linewidth',2);   
%    xP = tv; yP = Vs(tv);
%       plot(xP,yP,'k-','linewidth',2)       
%   
%    xlim([0 500]) 
%    ax = gca;ax.FontSize = FS;   
%    
%    plotMonths(xtxt, ytxt)
%    
%    hL = legend('I_{ds}','I_s','H_s','H_{ds}','V_s');
%    set(hL,'orientation','horizontal','location','northoutside','box','off')
%    
%    
%   
%   
% subplot('Position',[h2 v3 width height]) 
%   
%    xtxt = 'days';
%    ytxt = 'hospitalizations';grid on; box on;
%   
%  %  Hp = zeros(nT,1);
%  %  Hp(tp) = sf(end).*I(tp); 
%    xP = tp./10; yP = Hp;
%    plot(xP,yP,'r-','linewidth',2)
%    hold on
%    plot(ts,Hd(ts),'b-','linewidth',2)
%   
%    xlim([0 500]) 
%    yMax = ylim;
%    ylim([yMax(1) yMax(2)])
%    grid on; box on;
%    ax = gca;ax.FontSize = FS;
%    
%    xlabel(xtxt,'fontname','times','Fontsize',FS);
%    ylabel(ytxt,'fontname','times','Fontsize',FS);
% 
%    plotMonths(xtxt, ytxt)
%   
%    
% %    yyaxis right
% %    xP = ts; yP = 100*sf;
% %    plot(xP,yP,'k-','linewidth',2)
% %    
% 
%    
% figure(8)  % 88888888888888888888888888888888888888888888888888888888888
%    set(gcf,'units','normalized');
%    set(gcf,'Position', [0.5 0.2 0.25 0.25])
%    set(gcf,'color','w');
%    FS = 12;
%    
%    xP = tDC; yP = DC;
%    plot(xP,yP,'r','linewidth',2)
%    hold on
%    temp = find(DCd == 0, 1)-1; 
%    xP = 1:temp; yP = DCd(xP);
%    plot(xP,yP,'b','linewidth',2)
%    grid off; box on
%    plotMonths(xtxt, ytxt)
%    
%    xtxt = 'days';
%    ytxt = 'New Daily Cases';
%    xlabel(xtxt,'fontname','times','Fontsize',FS);
%    ylabel(ytxt,'fontname','times','Fontsize',FS);
%   

figure(11)  % 11 11 11 11 11 11 11 11 11 11 11 11 11 11 
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.05 0.05 0.35 0.8500])
   set(gcf,'color','w');
   FS = 12;
   width = 0.35; height = 0.20;
   v3 = 0.15; v2 = 0.45; v1 = 0.75;
   h1 = 0.10; h2 = 0.56;
   
% 1  I   
subplot('Position',[h1 v1 width height])
  xtxt = 'days';
  ytxt = 'I';
  xP = dayR; yP = Id; 
  plot(xP,yP,'b+')
  hold on
  xP = t; yP = I; 
  plot(xP,yP,'r','linewidth',2)
 
  txt = sprintf('   I_{tot} / pop = %3.1f %% \n',IP);
  xlim([0 500])
  yMax = ylim;
  ylim([0.6*yMax(1) yMax(2)])
  ax = gca;ax.FontSize = FS;
  
  Htitle = title(strcat(cn,{'     '},datestr(tC)));
  set(Htitle,'fontweight','normal','fontsize',10)
  plotMonths(xtxt, ytxt)
  
  
% 2 Itot  
subplot('Position',[h2 v1 width height])
   xtxt = 'days';
  ytxt = 'I_{tot}';
  xP = dayR; yP = Idtot; 
  plot(xP,yP,'b+')
  hold on
  xP = t; yP = Itot; 
  plot(xP,yP,'r','linewidth',2)
 
  xlim([0 500])
  yMax = ylim;
  ylim([0.6*yMax(1) yMax(2)])
  
  
  ylim([0 yMax(2)])
  
  ax = gca;ax.FontSize = FS;

  plotMonths(xtxt, ytxt)
  
  title(txt,'fontweight','norma','fontsize',10)
   
% 3 New
subplot('Position',[h1 v2 width height])
  xP = tDC; yP = DC;
   plot(xP,yP,'r','linewidth',2)
   hold on
   temp = find(DCd == 0, 1)-1; 
   xP = 1:temp; yP = DCd(xP);
   plot(xP,yP,'b','linewidth',1)
   grid off; box on

   yMax = ylim;
   ylim([0 yMax(2)])
   plotMonths(xtxt, ytxt)
   
   xtxt = 'days';
   ytxt = 'New Daily Cases';
   xlabel(xtxt,'fontname','times','Fontsize',FS);
   ylabel(ytxt,'fontname','times','Fontsize',FS);

% 4 D  
subplot('Position',[h2 v2 width height]) 
  xtxt = 'days';
  ytxt = 'D / 1000';
  xP = dayR; yP = Dd./1000; 
  plot(xP,yP,'b+')
  hold on
%   xP = t; yP = D./1000; 
%   plot(xP,yP,'r','linewidth',2)

   xP = t(10*Ndays:end); yP = D7(10*Ndays:end)./1000;
   plot(xP,yP,'r','linewidth',2)
   
  xlim([0 500])
  yMax = ylim;
  ylim([0.2*yMax(2) yMax(2)])
  
 % txt = sprintf('  D_{data} = %1.1f %%   D_{model} = %1.1f %% \n',Ddp, Dp);
   txt = sprintf('     ( D / I_{tot} )_{final} = %1.1f %% \n',100*Dindex);
%  title(txt,'fontweight','normal','fontsize',10)
  text(10,1.01*yMax(2),txt)
  ax = gca;ax.FontSize = FS;

  plotMonths(xtxt, ytxt)

% 5 RE  
subplot('Position',[h1 v3 width height]) 
   xtxt = 'days';
  ytxt = 'R_E';
  xP = t; yP = RE; 
  plot(xP,yP,'r','linewidth',2)
  hold on
  plot([0 500],[1 1],'k','linewidth',1)
  
  xlim([0 500])
  yMax = ylim;
  ylim([0 yMax(2)])
  
  ax = gca;ax.FontSize = FS;
  
  Htitle = title(strcat(cn,{'     '},datestr(tC)));
  set(Htitle,'fontweight','normal','fontsize',10)
  plotMonths(xtxt, ytxt)
  
% 6 DELTA  
subplot('Position',[h2 v3 width height]) 
   xtxt = 'days';
  ytxt = '\Delta';
  xP = t; yP = delta; 
  plot(xP,yP,'r','linewidth',2)
 
  xlim([0 500])
  yMax = ylim;
  ylim([0 yMax(2)])
  
  ax = gca;ax.FontSize = FS;
  plotMonths(xtxt, ytxt)


figure(9)  % 9999999999999999999999999999999999999999999999999999999999
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.5 0.2 0.35 0.35])
   set(gcf,'color','w');
   FS = 14;
   
    xtxt = 'days';
    ytxt = 'scaled pop.';grid on; box on;
   
   hold on
   xP = ts; yP = Ids(ts);
      plot(xP,yP,'g-','linewidth',2)
   xP = tp./10; yP = Isp;   
      plot(xP,yP,'r-','linewidth',2)
   xP = ts; yP = Hs(ts);
      plot(xP, yP,'b-','linewidth',2); 
   xP = tp./10; yP = Hps;   
      plot(xP,yP,'m-','linewidth',2)
    
   xP = tv; yP = Vs(tv);
      plot(xP,yP,'k-','linewidth',2)       
  
  % xlim([0 500]) 
   ax = gca;ax.FontSize = FS;   
   
   plotMonths(xtxt, ytxt)
   
   xlim([0 500]) 
   ylim([0 1.0]) 
   hL = legend('I_{dS}','I_s','H_{dS}','H_S','V_S');
   set(hL,'orientation','horizontal','location','northoutside','box','off')
   
   
  
toc

% =======================================================================

function  plotMonths(xtxt,ytxt) %plotMonths(z)
         hold on
         zz = [31 59 90 121 152 182 213 244 275 306 336 367 398 426 457 488];
         
        
         yMax = ylim;
         y2 = 1*yMax(2);
        
        y = get(gca,'ytick');
        y1 = 0;
        y3 = y(length(y));
        
        if y3 > y2; y2 = y3; end
        
        for c = 1 : length(zz)
        %  plot([zz(c) zz(c)],[0,z],'-','color',[0.5 0.5 0.5],'linewidth',0.01)
          plot([zz(c) zz(c)],[y1, 1.02*y2],'-','color',[0.5 0.5 0.5],'linewidth',0.01)
        end 
        
        Z2 = 'JFMAMJJASONDJFMA';
        
        
        for c = 1:16
          text(zz(c)-25,0.95*y2,Z2(c),'fontsize',8)
         %  text(zz(c)-25,z,Z2(c),'fontsize',8)
        end

        FS = 16;
        xlabel(xtxt,'fontname','times','Fontsize',FS);
        ylabel(ytxt,'fontname','times','Fontsize',FS);
       
        
        box on
        ax = gca;
        ax.XGrid = 'off';
        ax.YGrid = 'on';


end
  

