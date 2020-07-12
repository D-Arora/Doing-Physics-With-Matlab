% qm_PEB.m

% A numerical model for the Photoelectric Effect.
% Photocurretnas a function of stopping voltage and frequency

% Ian Cooper
% School of Physics, University of Sydney
% www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Notes:  http://www.physics.usyd.edu.au/teach_res/hsp/sp/mod7/m7PE.pdf


% 181028  Matlab 2018b

close all
clear 
clc

 
% SETUP ===============================================================
   lambdaR = [100, 800].*1e-9;      % wavelength limits   [m]    
   num = 500;
   e  = 1.60217662e-19;      % elementary charge  [C]
   me = 9.10938356e-31;      % electron mass      [kg]
   h  = 6.62607004e-34;      % Planck constant    [J.s]
   c  = 3.000e8;             % speed of light     [m/s]
   %Wmin = 4.70;              % Work Function      [eV]



% PHOTOCURRENT =========================================================
   nE = 10000;
   BE1 = 2; BE2 = 10*BE1;
   BE = (BE2 - BE1).*rand(nE,1) + BE1;
   fC = e*BE1/h;
   
   f = 5*fC;
   
   E = h*f/e;
   Int0 = [1,2,4];
   N = 200;
   vSmax = 10;
   vS = linspace(0,vSmax,N);
   I = zeros(20,1);
   Eesc = BE(BE<E) - BE1;
   Nanode = zeros(N,1);
   for c1 = 1:N
     Eanode = Eesc(Eesc > vS(c1));
     Nanode(c1) = length(Eanode);    
   end
   
   Nanode = Nanode ./ 1e3; 
  
figure(1)
  FS = 12;
  pos = [0.05 0.05 0.30 0.40];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
  
  xP = vS; yP = Int0(1).* Nanode;
    plot(xP,yP,'b','linewidth',2)
  hold on
  xP = vS; yP = Int0(2).* Nanode;
    plot(xP,yP,'r','linewidth',2)
  xP = vS; yP = Int0(3).* Nanode;
    plot(xP,yP,'m','linewidth',2) 
    
 ylim([0 20])   
  
  legend('1','2','4','orientation','horizontal') 
  text(2.1,19,'Incident Light Intensty','fontsize',FS)
  ylabel('current  I  [ a.u. ]')
  xlabel('stopping voltage  v_S  [ V ]')
  
   tm1 = '  f_C = ';
   tm2 = num2str(fC,'%2.2e\n');
   tm3 = '  Hz';
   tm4 = '      f / f_C = ';
   tm5 = num2str(f/fC,'%2.1f\n');
   tm = [tm1 tm2 tm3 tm4 tm5];
   h_title = title(tm,'fontsize',FS,'fontweight','normal');
  
  
  grid on
  box on
  set(gca,'fontsize',FS)


