% oscC00ND.m
% 190605  Matlab 2018b

% SIMULATION OF THE MOTION OF N COUPLED OSCILLATORS
% Lattice of N atoms is exciting by driving the displacement of atom #1.
% All parameters are changed within the Script
% [default values and units]

% DOING PHYSICS ONLINE 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Scripts
%     http://www.physics.usyd.edu.au/teach_res/mp/mscripts/
% Documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/oscC00NA.htm
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney


clear
close all
clc
tic

% ANIMATION and GIF file ==============================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flagAG
    flagAG = 1;    
% file name for animated gif   
    ag_name = 'ag_osc.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.2;
    frame1 = 0;

    
% SETUP =============================================================== 

% INPUT: Driving Force  frequency/ amplitude / excited atom 
   fD = 1.0;
   AD = 10;
   nD = 2;
   
% INPUT: Atoms for Fourier Transform and plots  [n1 = 5 n2 = 8]
   n1 = 5; n2 = 8;
 
% INPUT: Atoms (number N = 71 and mass m1 = 0.1  m2 = 0.2)
%        separation distance between atoms d = 1 /Spring constant kS = 10 
   N = 61;
   m1 = 0.1;  m2 = 0.2;
   d = 1;
   kS0 = 10;
      
% INPUT: Time steps Nt = 251 (odd number) / Max time interval tMax = 5
   Nt = 351;   % must be an ODD number   Nt = 251
   tMax = 6;
  
% INPUT: Initial displacement of atoms
   e = zeros(N,Nt);
  
% INPUT:  Damping constant  (zero dramping force b = 0)
   b = 0;
  
% Initial conditions    
%   e(1,1) = 0.10;
%   e(1,2) = 0.10;
 %  e(3,1) = -0.0025;
 %  e(3,2) = -0.0025;
    

% CALCULATIONS ======================================================== \

%  Masses and spring constants
   m = m1 .* ones(N,1); 
   kS = kS0 .* ones(N-1,1);
   m(2:2:N-1) = m2;
   m1 = m(1);
   m2 = m(2);
   
% Time span t  / time step dt
   t = linspace(0,tMax,Nt);
   dt = t(2) - t(1); 
   
% Driving force
  FD = zeros(N,Nt);
  FD(nD,:) =  AD*sin(2*pi*fD*t);
   
% Angular frequency w (omega) 
   wD = 2*pi*fD;
  

% Max phase velocity v0
%   v0 = d*sqrt(kS0/m0);
% Max frequencies omega0 (w0)  f0      
%   w0 = sqrt(4*kS0/m0);
%   f0 = w0 / (2*pi);
   
% Calculation constants
  KS = dt^2.*kS;   Kb = b*dt;   KD = dt^2;    
  
% FINITE DIFFERENCE METHOD: Displacements from equilibrium ==============
 
for nt = 1 : Nt-2
  for c = 2 : N-1
     e(c,nt+2) = 2*e(c,nt+1) - e(c,nt)...
                -(KS(c-1)/m(c)) * (e(c,nt+1) - e(c-1,nt+1)) ...
                -(KS(c)/m(c))   * (e(c,nt+1) - e(c+1,nt+1)) ...   
                - (Kb/m(c))     * (e(c,nt+1) - e(c,nt))     ...
                + (KD/m(c))     * FD(c,nt+1);     
  end
end

% Range for displacements   
eMax = zeros(N,1);
eMin = zeros(N,1);
for c = 1 : N
   eMax(c) = max(e(c,:));
   eMin(c) = min(e(c,:));
end   
   

% FOURIER TRANSFORMS ==================================================

% Frequency estimates fPeaks1 and fPeaks2
  fMin = 0.01;
  fMax = 2*fD;
  Nf = 2901;
  H1 = zeros(1,Nf); H2 = H1;
  f = linspace(fMin,fMax,Nf);
  nE = ceil(N/2);
  
  for c = 1:Nf
    g = e(n1,:) .*  exp(1i*2*pi*f(c)*t);
    H1(c) = simpson1d(g,0,tMax);
  end
  
  for c = 1:Nf
    g = e(n2,:) .*  exp(1i*2*pi*f(c)*t);
    H2(c) = simpson1d(g,0,tMax);
  end
    
  PSD1 = 2.*conj(H1).*H1;
  PSD1 = PSD1./max(PSD1);
  [~, yy] = findpeaks(PSD1,f,'MinPeakProminence',0.6);
  fPeaks1 = yy;
  
  PSD2 = 2.*conj(H2).*H2;
  PSD2 = PSD2./max(PSD2);
  [~, yy] = findpeaks(PSD2,f,'MinPeakProminence',0.6);
  fPeaks2 = yy;
  
% Wavelength estimate lambda
  clear g
  H3 = zeros(1,Nf); lambda = linspace(2,2*10/fD,Nf);
  xf = 0:N-1;
  fn = e(1:N,Nt)';
  for c = 1 : Nf
    g = fn .*  exp(1i*2*pi*xf/lambda(c));
    H3(c) = simpson1d(g,0,N-1);
  end
  PSD3 = 2.*conj(H3).*H3;
  PSD3 = PSD3./max(PSD3);
  [~, yy] = findpeaks(PSD3,lambda,'MinPeakProminence',0.6);
  wL_Peaks = yy;

 
  
% DISPERSION RELATIONSHIPS ============================================ 
% omega Plus wP  /  omega Minus wM
Nk = 501;
k = linspace(0,2*pi/(2*d),Nk);
wP = zeros(Nk,1); wM = wP;
for ck = 1: Nk
  wP(ck) = sqrt((kS0*(m1+m2)/(m1*m2)) * (1 + sqrt(1-4*m1*m2*sin(k(ck)*d)^2/(m1+m2)^2)));
  wM(ck) = sqrt((kS0*(m1+m2)/(m1*m2)) * (1 - sqrt(1-4*m1*m2*sin(k(ck)*d)^2/(m1+m2)^2)));
end
 wP_max = sqrt(2*kS0*(m1+m2)/(m1*m2));
 wP_min = sqrt(2*kS0/m1);
 wM_max = sqrt(2*kS0/m2);
 wM_min = 0;
 
 fP = wP./(2*pi);
 fP_min = min(fP);
 fP_max = max(fP);
 fM = wM./(2*pi);
 fM_min = min(fM);
 fM_max = max(fM);

% Angular frequency (omega_Peaks) of frequency peak 
  w_Peaks = 2*pi*fPeaks1;
 
% Estimate propagation constant k_Peak from fP
%   xx = find(fP < f1,1);
   k_Peaks = 2*pi / wL_Peaks(1); 
   
% Propagation velocity from spectrum values v_Peak = f(spectrum) * wL(spectrum)
   v_Peaks = fPeaks1*wL_Peaks(1);   

 
% GRAPHICS ============================================================

%% ANIMATION ===========================================================    
if flagAG == 1  
    % May have to adjust limits for different simulations
      Xmax = N*d; Xtick = 0:20:Xmax;     
    %  Xmax = 50; Xtick = 0:10:Xmax; 
    
% Lattice Displacements 9999999999999999999999999999999999999999999999
figure(9)
    pos = [0.05 0.5 0.4 0.4];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
    xlim([0 Xmax+1])
 %   set(gca,'xtick',Xtick)

for c = 1: 5 : Nt
  subplot('position',[0.1 0.8 0.8 0.1])   
    for cm = 1 : N
       xP = (cm*d + e(cm,c)); yP = 0;
       hPlot = plot(xP,yP,'ob');
       set(hPlot,'markersize',5,'markerfacecolor','b')
       hold on
    end
    
     xP = n1 + e(n1,c); yP = 0;
       hPlot = plot(xP,yP,'or');
       set(hPlot,'markersize',8,'markerfacecolor','r')
       
     xP = n2 + e(n2,c); yP = 0;
       hPlot = plot(xP,yP,'ok');
       set(hPlot,'markersize',8,'markerfacecolor','k')
       
    xP = nD + e(nD,c); yP = 0;
       hPlot = plot(xP,yP,'om');
       set(hPlot,'markersize',8,'markerfacecolor','m')   
       
       txt = ['t  =  ' num2str(t(c), '%2.2f')];
       title(txt)
       xlabel('x  [mm]')
       xlim([0 Xmax+1])
       
    %  set(gca,'xtick',Xtick)
      ylim([-1 1])
      grid on
      set(gca,'ytick',[])
      set(gca,'fontsize',12) 
      pause(0.01)
      hold off
      xlim([0 Xmax+1])
  %    set(gca,'xtick',Xtick)
      
 subplot('position',[0.1 0.2 0.8 0.4])  
   
    for cm = 1 : N
       xP = cm*d; yP = e(cm,c);
       hPlot = plot(xP,yP,'ob');
       set(hPlot,'markersize',6,'markerfacecolor','b')
       hold on
       xP = 1:N; yP = e(:,c);
       plot(xP,yP,'b','linewidth',2);  
       
       xP = n1; yP = e(n1,c);
       hPlot = plot(xP,yP,'or');
       set(hPlot,'markersize',6,'markerfacecolor','r')
       
       xP = n2; yP = e(n2,c);
       hPlot = plot(xP,yP,'ok');
       set(hPlot,'markersize',6,'markerfacecolor','k')
       
       xP = nD; yP = e(nD,c);
       hPlot = plot(xP,yP,'om');
       set(hPlot,'markersize',6,'markerfacecolor','m')
       
       grid on;
    end
      xlabel('x  [mm]')
      ylabel('left  e   right')
     
      xlim([0 Xmax+1])
   %   set(gca,'xtick',Xtick)
      ylim([-10 10])
      grid on
      %set(gca,'ytick',[])
      ylim([-1.2*max(eMax) 1.2*max(eMax)])
      set(gca,'fontsize',12) 
      grid on;
      
      pause(0.01)
      
      hold off
      xlim([0 Xmax+1])
    %  set(gca,'xtick',Xtick)
         
         frame1 = frame1 + 1;
         frame = getframe(9);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
         % On the first loop, create the file. In subsequent loops, append.
         if frame1 == 1
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end


end

end


%% Time evoltion of displacement of 2 atoms 11111111111111111111111111111
figure(1)
    pos = [0.02 0.1 0.3 0.25];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
   
    xP = t; yP =  e(2,:);
    plot(xP,yP,'b','linewidth',2)
    hold on
    yP =  e(n1,:);
    plot(xP,yP,'r','linewidth',2)
    yP =  e(n2,:);
    plot(xP,yP,'k','linewidth',2)
    legend(num2str(nD,'%2.0f'),num2str(n1,'%2.0f'),num2str(n2,'%2.0f'),'location','south','orientation','horizontal')
    xlabel('t  [ s ]')
    ylabel('left   e   right')
    grid on
    set(gca,'fontsize',12)
    ylim([-2.1 2.1])
    tm = sprintf('f_D = %3.2f   \n',fD);
    title(tm)
    

%% Fourier Transforms 2222222222222222222222222222222222222222222222222
figure(2)   
    pos = [0.33 0.1 0.3 0.5];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on  
subplot(2,1,1)   % frequency
    xP = f; yP = PSD1;
    plot(xP,yP,'r','linewidth',2)
    hold on
    yP = PSD2;
    plot(xP,yP,'k','linewidth',2)
    xlabel('f ')
    ylabel('PSD')
    grid on
    set(gca,'fontsize',12) 
    tm = sprintf('f_{Peak1} = %3.2f     f_{Peak2} = %3.2f \n',max(fPeaks1),max(fPeaks2));
    title(tm)
    set(gca,'fontsize',12)
    legend(num2str(n1,'%2.0f'),num2str(n2,'%2.0f'),'location','east')
    
subplot(2,1,2)    % wavelength
    xP = lambda; yP = PSD3;
    plot(xP,yP,'b','linewidth',2)
    xlabel('\lambda  ')
    ylabel('PSD [ a.u. ]')
    tm = sprintf('\\lambda_{Peak} = %3.2f \n',max(wL_Peaks));
    title(tm)
    grid on
    set(gca,'fontsize',12)
    

%%
figure(5)
    pos = [0.02 0.5 0.3 0.38];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
   
    xP = k.*(2*d/pi); yP =  wP./(2*pi);
    plot(xP,yP,'b','linewidth',2)
    hold on
    
    yP = wP_min.*ones(Nk,1)./(2*pi);
    plot(xP,yP,'b')
    
    yP =  wM./(2*pi);
    plot(xP,yP,'r','linewidth',2)
    
    yP = wM_max.*ones(Nk,1)./(2*pi);
    plot(xP,yP,'r')
    
    legend('f_+','f_+(min)','f_-','f_-(max)','location','south','orientation','horizontal')
    
    
    
    xlabel('k  [ \pi/2d ]')
    ylabel('f')
    grid on
    set(gca,'fontsize',12)
   
   tm1 = sprintf('f_-(min) = 0');
   tm2 = sprintf('  f_-(max) = %3.2f',wM_max./(2*pi));
   tm3 = sprintf('  f_+(min) = %3.2f',wP_min./(2*pi));
   tm4 = sprintf('  f_+(max) = %3.2f',wP_max./(2*pi));
   tm = [tm1 tm2 tm3 tm4]; 
   title(tm)
    
    
    %% Numerical output ==================================================== 
figure(4)   
    pos = [0.65 0.1 0.3 0.5];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on    
    
   xlim([0 100])
   ylim([30 160])
   fs = 14;
   
  
   Space = 160; space = 10;
     txt = ['Model Parameters:   atoms N = ' num2str(N,'%2.0f')];
     Htext = text(2,Space,txt);
     set(Htext,'fontsize',fs,'color','k')
  
   Space = Space - space; 
      txt = ['  Time steps N_t = ' num2str(Nt,'%3.0f')   '   Simulation time t_{Max} = ' num2str(tMax,'%3.2f')   ];
      Htext = text(2,Space,txt);
      set(Htext,'fontsize',fs,'color','k')
      
   Space = Space - space; 
      txt = ['   d = ' num2str(d,'%3.2f') '    m_1 = ' num2str(m1,'%3.2f') '    m_2 = ' num2str(m2,'%3.2f') '   k_{S0} = ' num2str(kS0,'%3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')
      
     Space = Space - space; 
      txt = ['   Driving force:  f_D = ' num2str(fD,'%3.2f') '   A_D = ' num2str(AD,'%3.2f')  ];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')   
      

  Space = Space - 2*space; 
      txt = ['Fourier Transforms - System response: Atom ' num2str(n1,'%2.0f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')   
     
 Space = Space - space; 
      txt = ['    f_{peak} = ' num2str(fPeaks1,'%3.2f ') ... 
         '    \omega_{peak} = ' num2str(w_Peaks,'%3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k') 

      Space = Space - space; 
      txt = ['   \lambda_{peak} = ' num2str(wL_Peaks,'%3.2f ') ... 
         '     k_{peak} = ' num2str(k_Peaks,'%3.2f') ' = ' num2str(k_Peaks*2*d/pi,'%3.2f') ' \pi/2d' ];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k') 
        
 Space = Space - space -2; 
      txt = ['    v_{peak} = ' num2str(v_Peaks,'%3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')    
      
% Space = Space - space - 5; 
%       txt = ['Theortical values given k = ' num2str(k_theory,'%3.2f') ...
%           ' = ' num2str(k_theory*d/pi,'%3.2f') '\pi/d' ];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')          
% 
%  Space = Space - space; 
%       txt = ['     \omega_{theory} = ' num2str(w_theory,'%2.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')     
%       
% Space = Space - space; 
%       txt = ['     f_{theory} = ' num2str(f_theory,'%2.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')        
%       
% Space = Space - space; 
%       txt = ['     vP_{theory} = ' num2str(vP_theory,'%2.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')         
%       
%  end
%  
% if df > 0
%  Space = Space - 2*space; 
%       txt = ['Fourier Transforms - System response: Atom ' num2str(n1,'%2.0f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')   
%      
%  Space = Space - space; 
%       txt = ['    f_{peak} = ' num2str(fPeaks1,'%3.2f    %3.2f \n') ... 
%              '   \lambda_{peak} = ' num2str(wL_Peaks,'%3.2f    %3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k') 
%    
% Space = Space - space; 
%       txt = ['    vP_1 = ' num2str(fPeaks1(1)*wL_Peaks(2),'%3.2f \n') ... 
%              '    vP_2 = ' num2str(fPeaks1(2)*wL_Peaks(1),'%3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k') 
%       
% end
% 
%       
  axis off   
%    

    toc

 