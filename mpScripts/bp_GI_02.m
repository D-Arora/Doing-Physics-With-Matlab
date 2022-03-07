% bp_GI_01.m


% Glucose - Insulin Dynamics

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/

% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/GImodelM.pdf

% IAN COOPER
% matlabvisualphysics@gmail.com

% 220211   Matlab R2021b

clear
%close all
clc

% Plot color
  col = [0 0 0]; 

% Model S M
  flagSM = 2;
   
  flagM = 2;

% TIME SCALE [minutes]  ================================================ 
  NT = 3600*24+1;
  tMax = 24*60;              % [1 day in minutes]
  t = linspace(0,tMax,NT);
  t = t';% [ minutes]
   % t = 60.*linspace(6,16,NT);
   % t = linspace(6*60,tMax+6*60,NT);  
  dt = t(2) - t(1);

% Change of units  mg/dL to mmol/L
%  s = 0.0555;


% MODEL PARAMETERS  =====================================================

% Basal production of glucose
   q10 = 2.1;
   q11 = 0.3;
   t11 = 3*60;
   t12 =1*60;

   dG_prod = q10.*ones(NT,1);
   if flagSM ==2
      dG_prod = dG_prod + q11.*exp( -(t - t11).^2 /(2.*t12^2) );
   end
 
% Glucose insulin independent glucose utilization
   q2 = 1.0e-3;
% Insulin sensitivity   
   q3 = 3e-3;


% Max rate of insulin secretion from pancreas
   Imax0 = 0.28;
   Imax = Imax0.*ones(NT,1);

   t21 = 15*60;
   t22 = 2*60;
   q21 =-0.04;
 
   if flagSM == 2
      Imax = Imax + q21.*exp( -(t - t21).^2 /(2.*t22^2) );
   end


% Sigmoidal function constant  
  a = 2e4; %1e4;
% Rate at which insulin is cleared (liver mainly)
  tHI = 69 ;%46;              % half-life  [min]
  kI = log(2)/tHI;
  kI = 0.01;

% Basal glucose and insulin (roots of cubic polynomial)
  p1 = q2+q3*Imax0/kI; p2 = -q10; p3 = q2*a; p4 = q3*Imax0/kI - q10*a;
  z = roots([p1 p2 p3 p4]);

  G_basal = z(1);
  I_basal = (Imax0/kI)*G_basal^2/(a+G_basal^2);


% FOOD INTAKE: STOMACH & GUT ==========================================

% GUT decay constant
    tH_GUT = 7;       % half-life   [min]
    k_GUT = log(2)/tH_GUT;

% MEALS (default is 6 meals)
  mL = 6;
% Quantity Q0
  G0 = [250 250 0 0 0 0];
% Time of meal [hours]  tM
  tM = [6 15 14 1 1 1];
% Meal half-life  [min]  
  tHG = [30 30 30 30 30 30]; 

  k = log(2)./tHG;

  G1 = zeros(NT,1);
  G2 = zeros(NT,1);
  G3 = zeros(NT,1);
  G4 = zeros(NT,1);
  G5 = zeros(NT,1);
  G6 = zeros(NT,1);

  R1 = tM(1)*3600:NT;
  R2 = tM(2)*3600:NT;
  R3 = tM(3)*3600:NT;
  R4 = tM(4)*3600:NT;
  R5 = tM(5)*3600:NT;
  R6 = tM(6)*3600:NT;

  G1(R1) = G0(1).*exp(-k(1) .* (t(R1)-t(R1(1))));
  G2(R2) = G0(2).*exp(-k(2) .* (t(R2)-t(R2(1))));     
  G3(R3) = G0(3).*exp(-k(3) .* (t(R3)-t(R3(1))));    
  G4(R4) = G0(4).*exp(-k(4) .* (t(R4)-t(R4(1))));   
  G5(R5) = G0(5).*exp(-k(5) .* (t(R5)-t(R5(1))));      
  G6(R6) = G0(6).*exp(-k(6) .* (t(R6)-t(R6(1))));      
    
    
% Stomach glucose    
    G_STO = G1 + G2 + G3 + G4 + G5  +G6;
    dG_STO = sum(G0);
    
% Rate of stomach glucose intake
    g_STO = k(1).*G1 + k(2).*G2 + k(3).*G3 + k(4).*G4 + k(5).*G5 + k(6).*G6 + k(6).*G6;


% Variation in glucose production over a day





% COMPUTATIONS  =======================================================

% Stomach glucose
  G_GUT = zeros(NT,1);
  dG_GUT = zeros(NT,1);
% Basal glucose production
   G_prod = zeros(NT,1);
  
% Glucose dclearance (decay);
   dG_decay = zeros(NT,1);
% insulin glucose clearance
   dG_insulin = zeros(NT,1);
% Blood glucose [mg/dL]  
  G = zeros(NT,1);
% Insulin  [uU/dL]
  I = zeros(NT,1);
  dI_pancrease = zeros(NT,1);
  dI_decay = zeros(NT,1);

  G(1) = G_basal;
  I(1) = I_basal;



for n = 1 : NT-1

% GLUCOSE -------------------------------------------------------------
% Rate of GUT glucose Intake
   G_GUT(n+1)  = G_GUT(n) + dt * (g_STO(n) - k_GUT*G_GUT(n));  
   dG_GUT(n) =  k_GUT*G_GUT(n);

%  Glucose production  
    if flagM == 2
        dG_prod(n) = 2*(dG_prod(n))/(1+exp(0.02*(G(n)-G_basal)));
    end

%  Glucose decay   
   dG_decay(n) = -q2*G(n);
   
%  Glucose metabolism by insulin
   dG_insulin(n) = - q3*I(n)*G(n);
   
%  Glucose total  
   G(n+1) = G(n) + dt*( dG_prod(n) + dG_insulin(n) + dG_decay(n) + dG_GUT(n) );


   % INSULIN  ------------------------------------------------------------    
% Insulin  pancrease
  dI_pancrease(n) =  Imax(n)*(G(n)^2/(a + G(n)^2));
  
% Insulin decay
  dI_decay(n) =  - kI*I(n); 

% INSULIN  
  I(n+1)    = I(n) + dt*( dI_pancrease(n) + dI_decay(n) );  
end

  dG_decay(NT) = dG_decay(NT-1);
  dG_insulin(NT) = dG_insulin(NT-1);
  dI_pancrease(NT) = dI_pancrease(NT-1);

%  
  G_avg = mean (G);
  [pks, locs] = findpeaks(G,'MinPeakHeight',1.1*G_basal);
  G_peaks = pks;
  G_locs  = t(locs)/60;
 
  % FWHM
%    G_half = G_basal + (G_peaks-G_basal)/2;
%    z1 = find(G>G_half,1);
%    z2 = find(G(locs:end) < G_half,1) + locs(1);
%    G_FWHM = (t(z2) - t(z1));
%  
  I_avg = mean(I);
  [pks, locs, width, z] = findpeaks(I,'MinPeakHeight',1.1*I_basal);
  I_peaks = pks;
  I_locs  = t(locs)/60;
 
%  FWHM  (comment / uncomment)
%    G_half = G_basal + (G_peaks-G_basal)/2;
%    z1 = find(G>G_half,1);
%    z2 = find(G(locs:end) < G_half,1) + locs(1);
%    G_FWHM = (t(z2) - t(z1))
% 
%    I_half = I_basal + (I_peaks-I_basal)/2;
%    z1 = find(I>I_half,1);
%    z2 = find(I(locs:end) < I_half,1) + locs(1);
%    I_FWHM = (t(z2) - t(z1))
   

G_load = (simpson1d(G',t(1), t(NT)) - G_basal*24*60)/(24*60);
I_load = (simpson1d(I',t(1), t(NT)) - I_basal*24*60)/(24*60);


% DISPLAY  ============================================================
disp(' ')
disp('Input Parameters')
fprintf('Basal production of glucose  a1 = %2.2f  mg/dL/min \n',q10)
fprintf('Glucose insulin independent glucose utilization  a2 = %2.2e  1/min \n',q2) 
fprintf('Insulin sensitivity   a3 = %2.3e  mL/uU/min \n',q3)
fprintf('Max rate of insulin secretion from pancreas   Imax = %2.2f mU/L/min \n',Imax0)
fprintf('Sigmoidal function constant   a = %2.3e  mg^2/dL^2  \n',a)
fprintf('Insulin clearance rate (liver mainly)half-life  tH_decay = %2.2f  min \n',tHI)
disp('  ')
fprintf('MEALS - Stomach: quantity = %2.2f \n',G0)
fprintf('MEALS - stomach: decay half-life = %2.2f  min \n',tHG)
fprintf('MEALS - stomach: total quantity = %2.2f \n',dG_STO)
fprintf('MEALS - gut: decay   half-life = %2.2f  min \n',tH_GUT)
disp('   ')
disp('Output Parameters')
fprintf('Glucose: basal level    G_basal =  %2.2f  mg/dL \n',G_basal)
fprintf('Glucose: average level  G_avg   =  %2.2f  mg/dL \n',G_avg)
fprintf('Glucose: peak levels    t_peaks =  %2.2f  h \n',G_locs)
fprintf('Glucose: peak levels    G_peak  =  %2.2f  mg/dL \n',G_peaks)
disp('  ')
% fprintf('Glucose: basal level    G_basal =  %2.2f  mmol/L \n',s*G_basal)
% fprintf('Glucose: average level  G_avg   =  %2.2f  mmol/L \n',s*G_avg)
% fprintf('Glucose: peak levels    G_peaks =  %2.2f  mmol/L \n',s*G_peaks)
%fprintf('Glucose: peak widths    G_FWHM =  %2.2f   min  \n',G_FWHM)
fprintf('Glucose: load    G_load =  %2.2f mg/dL  \n',G_load)
disp('  ')
fprintf('Insulin: basal level    I_basal =  %2.2f  mU/L \n',I_basal)
fprintf('Insulin: average level  I_avg   =  %2.2f  mU/L  \n',I_avg)
fprintf('Insulin: peak levels    t_peaks =  %2.2f  h \n',I_locs)
fprintf('Insulin: peak levels    I_peaks =  %2.2f  mU/L  \n',I_peaks)
%fprintf('Insulin: peak widths    I_FWHM =  %2.2f   min  \n',I_FWHM)
fprintf('Insulin: load    I_load =  %2.2f uU/mL  \n',I_load)


% GRAPHICS  ===========================================================

figure(9)  %  999999999999999999999999999999999999999999999999999999999
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.15 0.35 0.7]);
  set(gcf,'color','w');
  FS = 14;

subplot(4,1,1)
  xP = t./60; yP = dG_GUT;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([0 24])
%  ylim([0 0.1])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  ylabel('dG/dt_{GUT}')
  set(gca,'fontsize',FS)

subplot(4,1,2)
  xP = t./60; yP = dG_prod;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([0 24])
  ylim([0 2.5])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  ylabel('dG/dt_{prod}')
  set(gca,'fontsize',FS)

subplot(4,1,3)
  xP = t./60; yP = dG_decay;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([0 24])
%  ylim([-0.01 0])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  %set(gca,'ytick',0:10:max(y))
  ylabel('dG/dt_{decay}')
  set(gca,'fontsize',FS)

subplot(4,1,4)
  xP = t./60; yP = dG_insulin;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([0 24])
%  ylim([-0.1 0])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  ylabel('dG/dt_{insulin}')
  xlabel('time  t  [h]')
  set(gca,'fontsize',FS)


 figure(1)  % 111111111111111111111111111111111111111111111111111111111
  set(gcf,'units','normalized');
  set(gcf,'position',[0.2 0.20 0.35 0.7]);
  set(gcf,'color','w');
  FS = 14;
  
  sf = 0;
 % sf = s;   % mg/dL to mmol/L
  s = 1;
  
  limX = [0 24];

subplot(4,1,1)
  hold on
  xP = t./60; yP = G_STO;
  plot(xP,yP,'color',col,'linewidth',2)
  grid on; box on
 xlim(limX)
  set(gca,'xtick',0:2:24)
 % set(gca,'ytick',0:50:200)
  grid(gca,'minor')
  ylabel('G_{STO}')
  set(gca,'fontsize',FS)
  tm = sprintf('Stomach glucose G_{STO} =  %3.0f',dG_STO);
 % title(tm,'fontweight','normal')
 % text(14,240,'t_{STO} = 30  (high GI)','color','r','fontsize',FS);
 %  text(14,160,'t_{STO} = 60  (low GI)','color','b','fontsize',FS);
  
subplot(4,1,2)
hold on
  xP = t./60; yP = G_GUT;
  plot(xP,yP,'color',col,'linewidth',2)
  grid on; box on
  xlim(limX)
  set(gca,'xtick',0:2:24)
  grid(gca,'minor')
  ylabel('G_{GUT}')
  set(gca,'fontsize',FS)

subplot(4,1,3)
hold on
  xP = t./60; yP = s*G;
  plot(xP,yP,'color',col,'linewidth',2)
  hold on
  plot([0 t(end)/60],[G_basal G_basal],'m')
%   plot([0 t(end)/60],s.*[125 125],'m')
  grid on; box on
  xlim(limX)
%  ylim(s.*[60 260])
  set(gca,'xtick',0:2:24)
%  set(gca,'ytick',sf.*(0:25:150))
  grid(gca,'minor')
  ylabel('Glucose [mg.dL^{-1}]')
  if sf > 1; ylabel('Glucose   [mmolL^{-1}]'); end
 % txt = sprintf('G_{baslal} = %2.1f      G_{avg} = %2.1f', s*G_basal, s*G_avg);
 % title(txt,'fontweight','normal')
  set(gca,'fontsize',FS)

subplot(4,1,4)
hold on
  xP = t./60; yP = I;
  yP = I_basal +  (I-I_basal).*16;
  plot(xP,yP,'color',col,'linewidth',2)
  hold on
   plot([0 t(end)/60],[I_basal I_basal],'m')
  grid on; box on
  xlim(limX)
%  ylim([5 15])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  ylabel('Insulin   [mU.L^{-1}]')
  xlabel('Time  t  [h]')
  txt = sprintf('I_{baslal} = %2.1f      I_{avg} = %2.1f', I_basal, I_avg);
 % title(txt,'fontweight','normal')
  set(gca,'fontsize',FS)


figure(2)  % 222222222222222222222222222222222222222222222222222222222222
  set(gcf,'units','normalized');
  set(gcf,'position',[0.6 0.60 0.25 0.25]);
  set(gcf,'color','w');
  FS = 14;
  
  xP = I; yP = G;
  plot(xP,yP,'b','linewidth',2)
  hold on
  Hplot = plot(I_basal,G_basal,'bo');
  set(Hplot, 'markersize',8,'markerfacecolor','g','markeredgecolor','g')
  grid on
 
  ylabel('G','color','m')
  xlabel('I','FontName','Times')

  txt = sprintf('G_{basal} = %2.1f    I_{basal} = %2.1f ', G_basal, I_basal);
  title(txt,'fontweight','normal')
  set(gca,'fontsize',FS)
 
% Nullclines
G_NC = linspace(60,200,501);
I_GNC = (q10 - q2.*G_NC)./(q3.*G_NC);

yP = G_NC; xP = I_GNC;
plot(xP,yP,'m','linewidth',2)

xP = Imax0.*(G_NC.^2./(a+G_NC.^2))./kI;
plot(xP,yP,'k','linewidth',2)


  % ================================================================
 

  f = linspace(0,200,2000);
  S = zeros(2000,1);
  for n = 1 : 2000
        S(n) = sigm(f(n),0,1e3);
  end   
 
  
% figure(2)
%    set(gcf,'units','normalized');
%   set(gcf,'position',[0.5 0.10 0.25 0.35]);
%   set(gcf,'color','w');
%   FS = 14;
%   
%   plot(f,S)
  

% FUNCTIONS  ==========================================================

function S = sigm(f,f1,f2)
 
  F = (f-f1)^2;
  S = F / (f2 + F);

end

