% qm_PE.m

% A numerical model for the Photoelectric Effect.


% Ian Cooper
% School of Physics, University of Sydney
% Documentation: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
%                
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts


% 181027  Matlab 2018b

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


% =====================================================================   
% Stopping voltage vs frequency
% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flagS = 1 (save)
   flagS = 0;
%  Enter file name
     ag_name = 'agQM1.gif';
%  Delay in seconds before displaying the next image  
    delay = 2;  
%  Frame counter start
    nt = 1;

  fR = c./lambdaR;              % frequency limits    [Hz]
  f = linspace(0,fR(1),num);    % frequency range      [Hz]
  numW = 15;
  Wmin = linspace(1,8,numW);      % Work function [eV]
  
figure(1)
  FS = 14;
  pos = [0.05 0.05 0.35 0.5];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 

for cc = 1 : numW  
  vS = (h/e).*f - Wmin(cc);         % Stopping Voltage     {V}
  fC = e*Wmin(cc)/h;
   
  xP = f; yP = vS; 
    plot(xP,yP,'r','linewidth',1)
  hold on
 
  xP = f(vS>0); yP = vS(vS>0);
    plot(xP,yP,'b','linewidth',2)
  
  xP = 0; yP = -Wmin(cc);
    Hplot = plot(xP,yP,'o');
    set(Hplot,'markersize',8,'markerFaceColor','r','markerEdgeColor','r')
  xP = f(find(vS>0,1)); yP = 0;
    Hplot = plot(xP,yP,'o');
    set(Hplot,'markersize',8,'markerFaceColor','b','markerEdgeColor','b')
  
  xlabel('frequency  f  [ Hz ]')
  ylabel('stopping voltage  v_S  [ V ]')
   tm1 = '  W_{min} = ';
   tm2 = num2str(Wmin(cc),'%2.2f\n');
   tm3 = '  eV';
   tm4 = '      f_C = ';
   tm5 = num2str(fC,'%2.2e\n');
   tm6 = '  Hz';
   tm = [tm1 tm2 tm3 tm4 tm5 tm6];
   h_title = title(tm,'fontsize',FS,'fontweight','normal');
  
  ylim([-8 8])
  xlim([0 3e15])
  grid on
  box on
  set(gca,'fontsize',FS)
     if flagS == 1
       frame = getframe(1);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
       nt = nt+1;
     end    
 % pause(0.5)
  hold off
end  


% =====================================================================
% PHOTOCURRENT
% =====================================================================
%%
nE = 1000;
BE1 = 2; BE2 = 10*BE1;
BE = (BE2 - BE1).*rand(nE,1) + BE1;
fC = e*BE1/h;
E = 2*h*fC/e;
nP = 20;
vS = linspace(0,E,20);
I = zeros(20,1);
%for c1 = 1:1
  %  if vS(c1) < E
    %  for c2 = 1:nP 
    
    Eesc = BE(BE<E)
    Eanode = Eesc(Eesc>3)
    length(Eanode)
    
    % I(c1) = I(c1) +  length(E-BE > vS(c1));
   %   end
 %   end
%end

figure(2)
plot(vS,I)




