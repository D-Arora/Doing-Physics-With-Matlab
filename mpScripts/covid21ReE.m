% covid21ReE.m
% MINIMAL MODEL FOR THE TIME EVOLUTION OF THE COVID-19 ENDEMIC
%  (1) AUSTRALIA    (2) UNITED KINGDOM

% AUSTRALIA DATA
% Download: https://www.worldometers.info/coronavirus/

% Download: https://covidlive.com.au/report/daily-vaccinations/aus
% Data Matrix  covid21M
%   col 1:  Itot
%   col 2:  Iactive
%   col 3:  Deaths
%   col 4:  Hospitalizations
%   col 5:  Vacinations
%   col 6:  Daily New Cases

% UK DATA
% Download: https://www.worldometers.info/coronavirus/
% Download: https://ourworldindata.org/covid-vaccinations

% Data Matrix  covid21M
%   col  7:  Itot
%   col  8:  Iactive
%   col  9:  Deaths
%   col 10:  Hospitalizations
%   col 11:  Vacinations
%   col 12:  Daily New Cases


% REFERENCES: 
%    Description      covidMinModel.docx   
%    EXCEL database   covid19.xlsx


% ====================================================================
clear 
close all
clc

tic

% Select COUNTRY:  AUSTRALIA 1  /  UK 2   =============================

  flagC = 1;


% LOAD DATA from saved file ============================================
    load('covid21M.mat')
    
   switch flagC 
    case 1
      cn = '      AUSTRALIA';
      pop = 25862611;           % Population 
      Ndays  = find(covid21M(:,1) == 0,1) - 1;   % Days 
      Idtot  = covid21M(1:Ndays,1);   % Total infections
      Id     = covid21M(1:Ndays,2);   % Current (active) infections
      Dd     = covid21M(1:Ndays,3);   % Deaths
      Hd     = covid21M(:,4);         % Hosipitalizations
      Vd     = covid21M(:,5);         % Number of vaccinations 
      DCd    = covid21M(:,6);         % Daily new cases of infection
       
     case 2  
      cn = '    U.K.';   
      pop = 68323609;          % Population 
      Ndays  = find(covid21M(:,7) == 0,1) - 1;   % Days 
      Idtot  = covid21M(1:Ndays,7);    % Total infections
      Id     = covid21M(1:Ndays,8);    % Current (active) infections
      Dd     = covid21M(1:Ndays,9);    % Deaths
      Hd     = covid21M(:,10);         % Hosipitalizations
      Vd     = covid21M(:,11);         % Number of vaccinations 
      DCd    = covid21M(:,12);         % Daily new cases of infection
   end

    Cd = Idtot - Id - Dd;    % ReCoveries
    Rd = Cd + Dd;            % Removals
    dayR = 1: Ndays;         % Days  1, 2, ... , Ndays  


% MODEL SETUP  ========================================================
  nT = 5000; tMax = 500;    % 5000 time steps / 500 days / 10 steps = 1 day
  t = linspace(0,tMax,nT);
  h = t(2)-t(1);            % time step

  I = zeros(nT,1);          % Active (current) infections
  R = zeros(nT,1);          % Removals
  RE = zeros(nT,1);         % Reproduction factor
  delta = zeros(nT,1);      % Reproduction factor INCREMENT
  a = zeros(nT,1);          % a parameter
  b = zeros(nT,1);          % b parameter

% Initail values for active infections and removals 
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
        
%  DEATHS INDEX 
      Dindex = zeros(nT,1);
      Ddindex = Dd./Idtot;
      k = 1;
      for c = 1:Ndays
            Dindex(k:k+9) = Ddindex(c);
            k = k + 10;     
      end
       Dindex(10*Ndays:end) = Dindex(k-1);


% Setting reproduction factor increment delta  ========================
switch flagC  
  
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
    delta(3000:3200) = 0.5e-2;
    delta(3200:3300) = 1e-2;
    delta(3300:3400) = 2e-2;
    delta(3400:3500) = 3e-2;
    delta(3500:3600) = 10e-2;
    delta(3600:3700) = 14e-2;
    delta(3700:3800) = 16e-2;
     
 
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
     delta(3400:3500) = 2e-2;
     delta(3500:3600) = 3.0e-2;
     delta(3600:3700) = 4.0e-2;
     delta(3700:3800) = 1.5e-2;

end

 
% *******************************************************************
% CALCULATIONS  
% *******************************************************************
   delta = smooth(delta,200);
for c = 1 : nT-1
     RE(c+1) = RE(c) - h*(a(c)*RE(c)*I(c)) + delta(c);
     I(c+1) = I(c) + h*(b(c)*(RE(c+1)-1)*I(c));
     R(c+1) = R(c) + h*b(c)*I(c);
end

% Total Infections
  Itot = I + R;
% Deaths  
  D = Dindex .* Itot; 
% ReCoveries      
  C = R - D;


% Percentage of population infected
  IP = 100*Itot(end)/pop;
% Dates: tS start date . tC last ddate for data 
  tS = datetime(2021,1,1);
  tC = tS + Ndays-1;
  
% Daily New Cases
   Its = smooth(Itot,7);
   tDC = 1:500;
   DC = zeros(500,1);
   for c = 2:500
     DC(c) = Its(10*c) - Its(10*((c-1)));  
   end
   DC(1) = DC(2);
   DC = smooth(DC,7);

% Scaled active infections
   maxI = max([max(Id), max(I)]);
   Ids = Id./maxI;    % Data: scaled current infections
   Is  =  I./maxI;    % Model: scaled curremt infections

%   DATA: HOSPITALIZATIONS
   
        
     ts1 = find(Hd == 0, 1) - 2; % Data:  time span
     if ts1 > Ndays; ts1 = Ndays; end
     ts = 1:ts1;         
    
%  Model prediction: Hospitalizations;
     tp1 = 10*(Ndays-2);         % Model: time span
     tp = tp1:nT; 
     Isp = I(tp)./maxI;
     sf = mean(Hd(ts1-7:ts1)./Id(ts1-7:ts1));
     Hp = I(tp).* sf;

     maxH = max([max(Hd), max(Hp)]);
     Hs  = Hd./maxH;    % Scaled data hospitalizations
     Hps = Hp./maxH;    % Scaled model predictions

%  Scaled Vaccinations    
    tv1 = find(Vd == 0, 1) - 1;
    tv = 1:tv1;
    Vs = Vd./max(Vd);   % Scaled vaccination population
    

%  ******************************************************************* 
%  GRAPHICS 
%  *******************************************************************

  Z1 = [31	59	90	121	152	182	213	244	275	306	336	367	398	426	457	488];
  Z2 = 'JFMAMJJASONDJFMA';


figure(11)  % 11 11 11 11 11 11 11 11 11 11 11 11 11 11 
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.05 0.05 0.45 0.8500])
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
  
   txt = sprintf('   (I_{tot} / pop)_{final} = %3.1f %% \n',IP);
  title(txt,'fontweight','normal','fontsize',FS)
   
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
   ytxt = 'I_{new}';
   xlabel(xtxt,'fontname','times','Fontsize',FS);
   ylabel(ytxt,'fontname','times','Fontsize',FS);
   set(gca,'fontsize',FS)

% 4 D  
subplot('Position',[h2 v2 width height]) 
  xtxt = 'days';
  ytxt = 'D / 1000';
  xP = dayR; yP = Dd./1000; 
  plot(xP,yP,'b+')
  hold on
%   xP = t; yP = D./1000; 
%   plot(xP,yP,'r','linewidth',2)

   xP = t(10*Ndays:end); yP = D(10*Ndays:end)./1000;
   plot(xP,yP,'r','linewidth',2)
   
 %  plot(t,D/1000,'m','linewidth',2)
  xlim([0 500])
  yMax = ylim;
  ylim([0.2*yMax(2) yMax(2)])
  
 % txt = sprintf('  D_{data} = %1.1f %%   D_{model} = %1.1f %% \n',Ddp, Dp);
 txt = sprintf('     ( D / I_{tot} )_{final} = %1.1f %% \n',100*D(end)/Itot(end));
 % title(txt,'fontweight','normal','fontsize',10)
 text(80,1.02*yMax(2),txt,'fontsize',FS)
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
  
%  Htitle = title(strcat(cn,{'     '},datestr(tC)));
%  set(Htitle,'fontweight','normal','fontsize',10)
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
   set(gcf,'Position', [0.52 0.05 0.35 0.40])
   set(gcf,'color','w');
   FS = 14;
   
    xtxt = 'days';
    ytxt = 'scaled pop.';grid on; box on;
   
   hold on
   xP = ts; yP = smooth(Ids(ts),10);
      plot(xP,yP,'r-','linewidth',2)

  xP = tv; yP = smooth(Vs(tv),10);
      plot(xP,yP,'k-','linewidth',2)  

   xP = ts; yP = smooth(Hs(ts),10);
      plot(xP, yP,'b-','linewidth',2); 
  
   xP = tp./10; yP = smooth(Hps,10);   
      plot(xP,yP,'m-','linewidth',2)
    
       
  
   ax = gca;ax.FontSize = FS;   
 
   xlim([0 500])
   ylim([0 1.1])

   plotMonths(xtxt, ytxt)
   
   hL = legend('I_{dS}','V_{dS}','H_{dS}','H_{pS}');
   set(hL,'orientation','horizontal','location','northoutside','box','off')

figure(1)  % 11111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.52 0.55 0.45 0.30])
   set(gcf,'color','w');
   FS = 14;   
   
subplot(1,2,1)
   xP = t; yP = a;
   plot(xP, yP,'b-','linewidth',2); 
   ylim([0 1.1*max(a)])
   xlim([0 500])
   xticks(0:100:500)
   grid on
   xtxt = 'days';
   ytxt = 'a';
   xlabel(xtxt,'fontname','times','Fontsize',FS);
   ylabel(ytxt,'fontname','times','Fontsize',FS);
   set(gca,'fontsize',FS)


subplot(1,2,2)
   xP = t; yP = b;
   plot(xP, yP,'b-','linewidth',2); 
   ylim([0 1.1*max(b)])
   xlim([0 500])
   xticks(0:100:500)
   grid on

   xtxt = 'days';
   ytxt = 'b';
   xlabel(xtxt,'fontname','times','Fontsize',FS);
   ylabel(ytxt,'fontname','times','Fontsize',FS);
   set(gca,'fontsize',FS)


figure(2)  % 222222222222222222222222222222222222222222222222222222222
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.52 0.50 0.35 0.30])
   set(gcf,'color','w');
   FS = 14;   
   
   xP = ts; yP = Hd(ts);
     plot(xP, yP,'b-','linewidth',2); 
   hold on
  xP = tp./10; yP = Hp;
     plot(xP, yP,'r-','linewidth',2); 
 %  ylim([0.9*max(a) 1.1*max(a)])
   xlim([0 500])
   xticks(0:100:500)
   grid on
   xtxt = 'days';
   ytxt = 'H';
   xlabel(xtxt,'fontname','times','Fontsize',FS);
   ylabel(ytxt,'fontname','times','Fontsize',FS);

   yMax = ylim;
   ylim([0 yMax(2)])
  
   plotMonths(xtxt, ytxt)

   set(gca,'fontsize',FS)

   

% =======================================================================  
  toc
% =======================================================================


% FUNCTIONS ============================================================
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
  

