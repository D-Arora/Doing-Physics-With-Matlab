% oscC00ND.m
% 190605  Matlab 2018b

% SIMULATION OF THE MOTION OF N COUPLED OSCILLATORS
% Lattice of N atoms is exciting by driving the displacement of atom #1.
% All parameters are changed within the Script
% [default values and units]

% DOING PHYSICS ONLINE 
%    ../mphome.htm
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
    flagAG = 0;    
% file name for animated gif   
    ag_name = 'ag_osc.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.2;
    frame1 = 0;

    
% SETUP =============================================================== 

% INPUT: frequency f = 1 min >= 0.5 / frequency increment df = 0
   f1 = 2.00;
   df = 0.0;
% INPUT: Atoms for Fourier Transform and plots
   n1 = 11; n2 = 41;

% INPUT: Amplitude A of waves 1 and 2   
   A1 = 0.6; A2 = 0.6;
   
% INPUT: Atoms (number N = 71 and mass m = 0.1)
%        separation distance between atoms d = 1 /Spring constant kS = 10 
   N = 101;
   m = 0.1;
   d = 1;
   kS = 10;
      
% INPUT: Time steps Nt = 251 (odd number) / Max time interval tMax = 5
   Nt = 351;   % must be an ODD number   Nt = 251
   tMax = 10;
  
% INPUT: initial displacement of atoms
   e = zeros(N,Nt);
  
 
% CALCULATIONS ======================================================== \

% Time span t  / time step dt
   t = linspace(0,tMax,Nt);
   dt = t(2) - t(1); 
   
% Angular frequency w (omega) of waves 1 and 2 
  if f1 < 0.5, f1 = 0.5; end
   w1 = 2*pi*f1;
   f2 = f1 + df;
   w2 = 2*pi*f2;

% Max phase velocity v0
   v0 = d*sqrt(kS/m);
% Max frequencies omega0 (w0)  f0      
   w0 = sqrt(4*kS/m);
   f0 = w0 / (2*pi);
   
% Calculation constants
  KS = dt^2.*kS/m;       

  
% FINITE DIFFERENCE METHOD: Displacements from equilibrium ==============
 
for nt = 1 : Nt-2
  for c = 2 : N-1
    e(1,nt+1) = A1*sin(w1*(nt+1)*dt) + A2*sin(w2 *(nt+1)*dt) ;
      
    e(c,nt+2) = 2*e(c,nt+1) - e(c,nt)...
                -KS * (e(c,nt+1) - e(c-1,nt+1)) ...
                -KS * (e(c,nt+1) - e(c+1,nt+1));                
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
  fMax = 5;
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
  H3 = zeros(1,Nf); lambda = linspace(2,3*10/f1,Nf);
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

  
  
% Propagation parameters: velocity = frequency x wavelength
  vMax = d*sqrt(kS/m);
  wMax = sqrt(4*kS/m);
  fMax = wMax / (2*pi);
  
  if length(fPeaks1) == 1   
    v1 = wL_Peaks * fPeaks1;
    k_theory = 2*pi/max(wL_Peaks);
    w_theory = wMax*sin(k_theory * d/2);
    f_theory = w_theory / (2*pi);
    vP_theory = w_theory / k_theory;
  %  vG_theory = vMax * cos(k_theory*d/2);
  end
  


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
       hPlot = plot(xP,yP,'sb');
       set(hPlot,'markersize',4,'markerfacecolor','b')
       hold on
    end
    
     xP = n1 + e(n1,c); yP = 0;
       hPlot = plot(xP,yP,'or');
       set(hPlot,'markersize',6,'markerfacecolor','r')
       
       xP = n2 + e(n2,c); yP = 0;
       hPlot = plot(xP,yP,'ok');
       set(hPlot,'markersize',6,'markerfacecolor','k')
       
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
   
    xP = t; yP =  e(1,:);
    plot(xP,yP,'b','linewidth',2)
    hold on
    yP =  e(n1,:);
    plot(xP,yP,'r','linewidth',2)
    yP =  e(n2,:);
    plot(xP,yP,'k','linewidth',2)
    legend('1',num2str(n1,'%2.0f'),num2str(n2,'%2.0f'),'location','south','orientation','horizontal')
    xlabel('t  [ s ]')
    ylabel('left   e   right')
    grid on
    set(gca,'fontsize',12)
    ylim([-2.1 2.1])
    tm = sprintf('f_1 = %3.2f    f_2 = %3.2f   \n',f1,f2);
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
   % title(tm)
    set(gca,'fontsize',12)
    legend(num2str(n1,'%2.0f'),num2str(n2,'%2.0f'),'location','east')
    
subplot(2,1,2)    % wavelength
    xP = lambda; yP = PSD3;
    plot(xP,yP,'b','linewidth',2)
    xlabel('\lambda  ')
    ylabel('PSD [ a.u. ]')
  %  tm = sprintf('\\lambda_{Peak} = %3.2f     v = %3.2f   \n',max(wL_Peaks),v1);
    grid on
   % title(tm)
    set(gca,'fontsize',12)
    
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
      txt = ['   d = ' num2str(d,'%3.2f') '    m = ' num2str(m,'%3.2f') '    k_S = ' num2str(kS,'%3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')
      
  Space = Space - space; 
      txt = ['   \omega_{max} = ' num2str(w0,'%3.2f') '    f_{max} = ' num2str(f0,'%3.2f') '    v_{max} = ' num2str(v0,'%3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')
      
  Space = Space - space; 
      txt = ['   Driving frequencies  f_1 = ' num2str(f1,'%3.2f') '    f_2 = ' num2str(f2,'%3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')    
 
  if df == 0 
  Space = Space - 2*space; 
      txt = ['Fourier Transforms - System response: Atom ' num2str(n1,'%2.0f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')   
     
 Space = Space - space; 
      txt = ['    f_{peak} = ' num2str(fPeaks1,'%3.2f    %3.2f') ... 
         '    \omega_{peak} = ' num2str(2.*pi.*fPeaks1,'%3.2f    %3.2f') ...
         '   \lambda_{peak} = ' num2str(wL_Peaks,'%3.2f    %3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k') 

 
 Space = Space - space -2; 
      txt = ['    k_{peak} = ' num2str(k_theory,'%3.2f') ...
          '       v_P = f_{peak} \lambda_{peak} = ' num2str(v1,'%3.2f    %3.2f')  ];
      
      
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')    
      
Space = Space - space - 5; 
      txt = ['Theortical values given k = ' num2str(k_theory,'%3.2f') ...
          ' = ' num2str(k_theory*d/pi,'%3.2f') '\pi/d' ];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')          

 Space = Space - space; 
      txt = ['     \omega_{theory} = ' num2str(w_theory,'%2.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')     
      
Space = Space - space; 
      txt = ['     f_{theory} = ' num2str(f_theory,'%2.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')        
      
Space = Space - space; 
      txt = ['     vP_{theory} = ' num2str(vP_theory,'%2.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')         
      
 end
 
if df > 0
 Space = Space - 2*space; 
      txt = ['Fourier Transforms - System response: Atom ' num2str(n1,'%2.0f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k')   
     
 Space = Space - space; 
      txt = ['    f_{peak} = ' num2str(fPeaks1,'%3.2f    %3.2f \n') ... 
             '   \lambda_{peak} = ' num2str(wL_Peaks,'%3.2f    %3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k') 
   
Space = Space - space; 
      txt = ['    vP_1 = ' num2str(fPeaks1(1)*wL_Peaks(2),'%3.2f \n') ... 
             '    vP_2 = ' num2str(fPeaks1(2)*wL_Peaks(1),'%3.2f')];
      Htext = text(2,Space,txt); 
      set(Htext,'fontsize',fs,'color','k') 
      
end

      
 axis off   
   

    toc

 