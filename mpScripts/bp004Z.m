% bp_GI_01.m

% Gluocse - Insulin Dynamics

% bp_
clear
close all
clc


% plot color
  col = [0 0 1];
  
% TIME SCALE [minutes] 
  NT = 3600*24+1;
  tMax = 24*60;              % [1 day in minutes]
  t = linspace(0,tMax,NT);   % [ minutes]
  dt = t(2) - t(1);

% Change of units  mg/dL to mmol/L
  s = 0.0555;
  
% MODEL PARAMETERS  =====================================================
%  f = 1/(24*60);
  
% Basal production of glucose   [2.1]
   a1 = 2.1;
 %  a1 = 864*f;
% Glucose insulin independent glucose utilization   [1.0e-3]
   a2 = 1.0e-4;
 %  a2 = 1.44*f;
   
% Insulin sensitivity     [3.0e-3]
   a3 = 3.0e-3;
  % a3 = 5e-4;
 %  a3 = 0.72*f/(1+0.01*100);
   
% Max rate of insulin secretion from pancreas   [0.28]
  Imax = 0.28; 
 % Imax = 10;
%  Imax = 600*43.2*f;
  
% Sigmoidal function constant    [1e4]
  a = 1e4;
  a = 2e4;
% Rate at which insulin is cleared (liver mainly)    [46]
  tHI = 46;              % half-life  [min]
  kI = log(2)/tHI;
  
 % kI = 432*f;
  
% Gut decay constant
    tH_GUT = 7;       % half-life   [min]
    k_GUT = log(2)/tH_GUT;  

% Basal glucose and insulin (roots of cubic polynomial)
  p1 = a2+a3*Imax/kI; p2 = -a1; p3 = a2*a; p4 = a3*Imax/kI - a1*a;
  z = roots([p1 p2 p3 p4]);

  G_basal = z(1);
  I_basal = (Imax/kI)*G_basal^2/(a+G_basal^2);

 
  
% FOOD INTAKE: STOMACH & GUT ==========================================
       
% Meal #1
    G0(1)   = 200;
    tHG(1)  = 30;       % t_STO half-life  [min]
    
    k(1) = log(2)/tHG(1);
    G1 = zeros(NT,1);
    R1 = 8*3600:NT;
    G1(R1) = G0(1).*exp(-k(1) .* (t(R1)-t(R1(1))));
    
    G_STO = G1;
    dG_STO = sum(G0);
    g_STO = k(1)*G1;
%  % Meal #2
%     G0(2) = 0;
%     tHG(2) = 20;     % half-life  [min]
%     
%     k(2) = log(2)/tHG(2);
%     G2 = zeros(NT,1);
%     R2 = 10*3600:NT;
%     G2(R2) = G0(2).*exp(-k(2) .* (t(R2)-t(R2(1))));

% Stomach glucose    
%     G_STO = G1 + G2;
%     dG_STO = sum(G0);
%     
% % Rate of stomach glucose intake
%     g_STO = k(1)*G1 + k(2).*G2;
    

% COMPUTATIONS  =======================================================

% Stomach glucose
  G_GUT = zeros(NT,1);
  dG_GUT = zeros(NT,1);
% Basal glucose production
%  G_prod = t.*a1;
   dG_prod = zeros(NT,1);
% Glucose dclearance (decay);
   dG_decay = zeros(NT,1);
% insulin glucose clearance
   dG_insulin = zeros(NT,1);
% Blood glucose [mg/dL]  
  G = zeros(NT,1);
% Insulin  [uU/dL]
  I = zeros(NT,1);

  G(1) = G_basal;
  I(1) = I_basal;



for n = 1 : NT-1

% GLUCOSE -------------------------------------------------------------
% Rate of GUT glucose Intake
 
   G_GUT(n+1)  = G_GUT(n) + dt * (g_STO(n) - k_GUT*G_GUT(n));  
     
% Glucose production  
   dG_prod(n+1) = dt*a1;
    
% Glucose decay   
   dG_decay(n+1) = -dt*a2*G(n);
    
% Glucose metabolism by insulin
   dG_insulin(n+1) = - dt*a3*I(n)*G(n);
   
% Glucose total and time increment  
   G(n+1) = G(n) + dt*(a1 - a2*G(n) - a3*I(n)*G(n) + k_GUT*G_GUT(n));
   
   dG_GUT(n+1) = dt*k_GUT*G_GUT(n);

  
% INSULIN  ------------------------------------------------------------    
% Insulin  pancrease
  I_pancrease = dt * Imax*(G(n)^2/(a + G(n)^2)) ;

% I_decay
  I_decay =  - dt*kI*I(n);
    
  I(n+1)    = I(n) + I_pancrease + I_decay;  
end

 dG_prod(1) = dG_prod(end);
 dG_decay(1) = dG_decay(end);
 dG_insulin(1) = dG_insulin(end);

%  % Average, peak, FWHM, load
%    G_avg = mean (G);
%    [pks, locs] = findpeaks(G,'MinPeakHeight',1.1*G_basal);
%    G_peaks = pks;
%    G_locs  = t(locs)/60;
%  
%    G_half = G_basal + (G_peaks-G_basal)/2;
%    z1 = find(G>G_half,1);
%    z2 = find(G(locs:end) < G_half,1) + locs(1);
%    G_FWHM = (t(z2) - t(z1));
%  
%    I_avg = mean(I);
%    [pks, locs, width, z] = findpeaks(I,'MinPeakHeight',1.1*I_basal);
%    I_peaks = pks;
%    I_locs  = t(locs)/60;
%  
%    I_half = I_basal + (I_peaks-I_basal)/2;
%    z1 = find(I>I_half,1);
%    z2 = find(I(locs:end) < I_half,1) + locs(1);
%    I_FWHM = (t(z2) - t(z1));
%   
%   G_load = (simpson1d(G',t(1), t(NT)) - G_basal*24*60)/(24*60);
%   I_load = (simpson1d(I',t(1), t(NT)) - I_basal*24*60)/(24*60);  
%  % G_load = (simpson1d(G',t(1), t(NT)));
%  % I_load = (simpson1d(I',t(1), t(NT)));  
% % Time to return to within 5% of G_basal and I_basal [min]
% %   assuming meal at 0800
%   z = find(G>G_basal*1.05);
%   tGreturn = t(z(end)) - 8*60;
%   z = find(I>I_basal*1.05);
%   tIreturn = t(z(end)) - 8*60;
%   
%  
% DISPLAY  ============================================================
% disp(' ')
% disp('Input Parameters')
% fprintf('Basal production of glucose  a1 = %2.2f  mg/dL/min \n',a1)
% fprintf('Glucose insulin independent glucose utilization  a2 = %2.2e  L/min \n',a2) 
% fprintf('Insulin sensitivity   a3 = %2.3e  mL/uU/min \n',a3)
% fprintf('Max rate of insulin secretion from pancreas   Imax = %2.2f uU/mL/min \n',Imax)
% fprintf('Sigmoidal function constant   a = %2.3e  mg^2/dL^2  \n',a)
% fprintf('Insulin clearance rate (liver mainly) half-life  tH_decay = %2.2f  min \n',tHI)
% fprintf('GUT: decay   half-life = %2.2f  min \n',tH_GUT)
% 
% disp('  ')
% fprintf('MEALS - Stomach: quantity = %2.2f \n',G0)
% fprintf('MEALS - stomach: decay half-life = %2.2f  min \n',tHG)
% 
% disp('   ')
% disp('GLUCOSE Output Parameters')
% fprintf(' basal level   G_basal   =  %2.2f  mg/dL \n',G_basal)
% fprintf(' average level  G_avg    =  %2.2f  mg/dL \n',G_avg)
% fprintf(' peak level    G_peak    =  %2.2f  mg/dL \n',G_peaks)
% fprintf(' G_peak - G_basal        =  %2.2f mg/dL \n',G_peaks - G_basal)
% fprintf(' G_load                  =  %2.2f  mg/dL \n',G_load)
% fprintf(' peak widths   G_FWHM    =  %2.2f  min  \n',G_FWHM)
% fprintf(' peak time   dt_Gpeak    =  %2.2f  min \n',60*(G_locs - 8))
% fprintf(' return time dt_Greturn  =  %3.2f min \n',tGreturn)
% 
% disp('  ')
% disp('INSULIN Output Parameters')
% fprintf(' basal level   I_basal   =  %2.2f  uU./mL \n',I_basal)
% fprintf(' average level   I_avg   =  %2.2f  uU/mL  \n',I_avg)
% fprintf(' peak level    I_peaks   =  %2.2f  uU/mL  \n',I_peaks)
% fprintf(' I_peak - I_basal        =  %2.2f  mg/dL \n',I_peaks - I_basal)
% fprintf(' load           I_load   =  %2.2f  uU/mL \n',I_load)
% fprintf(' peak widths    I_FWHM   =  %2.2f  min  \n',I_FWHM)
% fprintf(' peak time   dt_Ipeak    =  %2.2f  min \n',60*(I_locs - 8))
% fprintf(' return time dt_Ireturn  =  %3.2f  min \n',tIreturn)
% disp('  ')
% %fprintf('Glucose: basal level    G_basal =  %2.2f  mmol/L \n',s*G_basal)
% %fprintf('Glucose: average level  G_avg   =  %2.2f  mmol/L \n',s*G_avg)
% %fprintf('Glucose: peak levels    G_peaks =  %2.2f  mmol/L \n',s*G_peaks)


% GRAPHICS  ===========================================================

 figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.10 0.35 0.7]);
  set(gcf,'color','w');
  FS = 13;
  
  sf = 1;
 % sf = s;   % mg/dL to mmol/L
%  box on
  
subplot(4,1,1)
hold on
  xP = t./60; yP = G_STO;
  plot(xP,yP,'color',[0 0 1],'linewidth',2)
  grid on
  xlim([0 24])
  set(gca,'xtick',0:2:24)
 % set(gca,'ytick',0:50:200)
  grid(gca,'minor')
  ylabel('G_{STO}')
  set(gca,'fontsize',FS)
  tm = sprintf('Stomach glucose G_{STO} =  %3.0f',dG_STO);
  box on
 % title(tm,'fontweight','normal')

subplot(4,1,2)
hold on
  xP = t./60; yP = G_GUT;
  plot(xP,yP,'color',[0 0 1],'linewidth',2)
  grid on
  xlim([0 24])
  set(gca,'xtick',0:2:24)
  grid(gca,'minor')
  ylabel('G_{GUT}')
  set(gca,'fontsize',FS)
  box on
  
subplot(4,1,3)
hold on
  xP = t./60; yP = sf*G;
  plot(xP,yP,'color',col,'linewidth',2)
  hold on
 % plot([0 t(end)/60],sf.*[70 70],'m')
 % plot([0 t(end)/60],sf.*[125 125],'m')
  grid on
  xlim([0 24])
 % xlim([6 16])
  ylim(sf.*[50 300])
  set(gca,'xtick',0:2:24)
%  set(gca,'ytick',sf.*(0:25:150))
  grid(gca,'minor')
  ylabel('Glucose [mg.dL^{-1}]')
  if sf > 1; ylabel('Glucose [mmolL^{-1}]'); end
%   txt = sprintf('G_{baslal} = %2.1f      G_{avg} = %2.1f', sf*G_basal, sf*G_avg);
%   title(txt,'fontweight','normal')
%    txt = sprintf('G_{basal} = %2.1f  mg.dL^{-1}',sf*G_basal);
%    text(15,155,txt,'fontsize',14)
  set(gca,'fontsize',FS)
  box on
  
subplot(4,1,4)
hold on
  xP = t./60; yP = I;
  plot(xP,yP,'color',col,'linewidth',2)
  grid on
  xlim([0 24])
 % xlim([6 16])
 % ylim([6 14])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  ylabel('Insulin [\muU.mL^{-1}]')
  xlabel('time  t  [h]')
 % txt = sprintf('I_{baslal} = %2.1f  \theta U.L_{-1}      I_{avg} = %2.1f', I_basal, I_avg);
 % title(txt,'fontweight','normal')
 % txt = sprintf('I_{basal} = %2.1f \\muU.L^{-1}  \n', I_basal);
 % text(15,10,txt,'fontsize',14)
 set(gca,'fontsize',FS)
 box on


figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.5 0.10 0.35 0.7]);
  set(gcf,'color','w');
  FS = 14;
  
subplot(4,1,1)
  xP = t./60; yP = dG_GUT;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([0 24])
  ylim([0 0.1])
  set(gca,'xtick',0:2:24)
 % set(gca,'ytick',0:50:200)
  grid(gca,'minor')
  ylabel('\DeltaG_{GUT}')
  set(gca,'fontsize',FS)
 

subplot(4,1,2)
  xP = t./60; yP = dG_prod;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([0 24])
  ylim([0 0.1])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  ylabel('\DeltaG_{prod}')
  set(gca,'fontsize',FS)

subplot(4,1,3)
  xP = t./60; yP = dG_decay;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([0 24])
  ylim([-0.1 0])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  %set(gca,'ytick',0:10:max(y))
  ylabel('\DeltaG_{decay}')
  set(gca,'fontsize',FS)


subplot(4,1,4)
  xP = t./60; yP = dG_insulin;
  plot(xP,yP,'b','linewidth',2)
  grid on
  xlim([0 24])
  ylim([-0.1 0])
  set(gca,'xtick',0:2:24)
   grid(gca,'minor')
  ylabel('\DeltaG_{insulin}')
  xlabel('time  t  [h]')
  set(gca,'fontsize',FS)



 

%   f = linspace(0,200,2000);
%   S = zeros(2000,1);
%   for n = 1 : 2000
%         S(n) = sigm(f(n),0,1e3);
%   end   
 
  
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

