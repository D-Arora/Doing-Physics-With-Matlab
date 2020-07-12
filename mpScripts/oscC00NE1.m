% oscC00NE.m
% 1900601  Matlab 2018b

% SIMULATION OF THE MOTION OF N COUPLED OSCILLATORS
% normal modes of vibration of an atomic lattice
% All parameters are changed within the Script
% [variables: default values and S.I. units (have large atoms)]
% Calls external function  simpson1d.m

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

% INPUT: [Atoms N = 61  / mass m = 0.1 kg  / spring constant kS 10 kg/m]
%        [separation distance between atoms d = 1 m] 
   N = 31;
   m0 = 0.1;
   kS0 = 10;
   d = 1;
   
% INPUT: [Time steps Nt = 351 s (odd number)/Max time interval tMax = 18 s]
   Nt = 551;   % must be an ODD number   
 %  tMax = 20;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
% Driving force parameters: Select an atom at an antinode for excitation
% INPUT: [frequency fD0 = 1/6 Hz  /  amplitude frequency AD = 10]
%         [nE atom exited by driving force]
%        [nE = n1 = 16 n2 = 11  atoms for Fourier Transform t/f calculations
  fD0 = 1/6;    
  AD = 0.08;
  nE = 16;      n1 = nE;
  n2 = 11;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    

% INPUT: [initial displacement of atoms e = 0]
   e = zeros(N,Nt);

% Damping constant   [b = 0.1  kg.s^-1] 
    b = 0.0;   
   

% CALCULATIONS ========================================================


% System parameters
  m = m0 .* ones(N,1); 
  kS = kS0 .* ones(N-1,1);

% Time span t  / time step dt
   tMax = 8/fD0;
   t = linspace(0,tMax,Nt);
   dt = t(2) - t(1); 

% Driving force    
   FD = zeros(N,Nt);
   FD(nE,:) =  AD*sin(2*pi*fD0*t); 

% Max phase velocity v0  [assume constant m and kS values]
   v0 = d*sqrt(kS0/m0);
% Max allowed propagation frequencies [omega0 (w0)  /  f0]      
   w0 = sqrt(4*kS0/m0);
   f0 = w0 / (2*pi);
   
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
  EMAX = max(max([abs(eMax), abs(eMin)])); 

% FOURIER TRANSFORMS ==================================================

% Frequency estimates for atoms n1 and n2 [fPeaks1 and fPeaks2]
%  Foureir function g  /  Fourier Transform H
%  INPUT: [Frequency range for Fouier Transform and number of grid points]
    fMin = 0.01; fMax = 3*fD0; Nf = 2901;
   
    f = linspace(fMin,fMax,Nf);
    H1 = zeros(1,Nf); H2 = H1;

% Calculate Fourier Transform
%         [power spectral density PSD / peaks in spectrum fPeaks]
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

% Calculate Fourier Transform for wavelength [lambda / peaks wL_Peaks]  
    clear g
    H3 = zeros(1,Nf); lambda = linspace(0.001,3*10/fD0,Nf);
    xf = 0:N-1;
    fn = e(1:N,Nt)';
  for c = 1 : Nf
    g = fn .*  exp(1i*2*pi*xf/lambda(c));
    H3(c) = simpson1d(g,0,40);
  end
    PSD3 = 2.*conj(H3).*H3;
    PSD3 = PSD3./max(PSD3);
    [~, yy] = findpeaks(PSD3,lambda,'MinPeakProminence',0.6);
    wL_Peaks = yy;

    
% OUTPUTS =============================================================
fprintf('Frequency peaks for n1  %3.4f \n',fPeaks1)
disp('   ')
fprintf('Frequency peaks for n2  %3.4f \n',fPeaks2)
disp('   ')
fprintf('Wavelegnth peaks  %3.4f \n',wL_Peaks)


% GRAPHICS ============================================================

% ANIMATION ===========================================================    
if flagAG == 1  
    % May have to adjust limits for different simulations
      Xmax = N*d; Xtick = 0:20:Xmax;     
    %  Xmax = 50; Xtick = 0:10:Xmax; 
    
% Lattice Displacements 9999999999999999999999999999999999999999999999
figure(9)
    pos = [0.05 0.5 0.3 0.4];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on
    xlim([0 Xmax+1])
 %   set(gca,'xtick',Xtick)

for c = 1: 5 : Nt
  subplot('position',[0.15 0.8 0.8 0.1])   
    for cm = 1 : N
       xP = (cm*d + e(cm,c)); yP = 0;
       hPlot = plot(xP,yP,'ob');
       set(hPlot,'markersize',8,'markerfacecolor','b')
       hold on
    end
       xP = [1 N]; yP = [0 0];
       hPlot = plot(xP,yP,'or');
       set(hPlot,'markersize',8,'markerfacecolor','r')
       
     
    %  xlabel('x  [m]')
      xlim([0 Xmax+1])
    %  set(gca,'xtick',Xtick)
      ylim([-1 1])
      grid on
      set(gca,'ytick',[])
      set(gca,'fontsize',12) 
      pause(0.01)
      hold off
      xlim([0 Xmax+1])
      %set(gca,'xtick',Xtick)
      
 subplot('position',[0.15 0.2 0.8 0.4])  
      ylim([-1.25 1.25])      

    for cm = 1 : N
       xlim([0 Xmax+1])
       ylim([-1.25 1.25])  
       xP = cm*d; yP = e(cm,c);
       hPlot = plot(xP,yP,'ob');
       set(hPlot,'markersize',8,'markerfacecolor','b')
       hold on
       ylim([-1.25 1.25])  
       xP = 1:N; yP = e(:,c);
       plot(xP,yP,'b','linewidth',2);  
       xP = [1 N]; yP = [0 0];
       hPlot = plot(xP,yP,'or');
       set(hPlot,'markersize',8,'markerfacecolor','r')
        xP = n1; yP = e(n1,c);
        hPlot = plot(xP,yP,'om');
        set(hPlot,'markersize',8,'markerfacecolor','m')
        xP = n2; yP = e(n2,c);
        hPlot = plot(xP,yP,'ok');
        set(hPlot,'markersize',8,'markerfacecolor','k') 
        
       grid on;
    end
      xlabel('x  [m]')
      ylabel('left  e   right')
      title(num2str(t(c),'%2.2f'))
      xlim([0 Xmax+1])
      ylim([-1.25 1.25])
      set(gca,'fontsize',12) 
      grid on;
      pause(0.0001)
      hold off
     
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


% Time evoltion of displacement of 2 atoms 11111111111111111111111111111
figure(1)
    pos = [0.02 0.1 0.3 0.25];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    box on
   
    xP = t; yP = e(n1,:);
    plot(xP,yP,'m','linewidth',2)
    hold on
    yP =  e(n2,:);
    plot(xP,yP,'k','linewidth',2)
    legend(num2str(n1,'%2.0f'),num2str(n2,'%2.0f'),'orientation','horizontal','location','northwest')
    xlabel('t  [ s ]')
    ylabel('left   e   right')
    grid on
    set(gca,'fontsize',12)
%    ylim([-2.1 2.1])
%    tm = sprintf('f_1 = %3.2f    f_2 = %3.2f   \n',f1,f2);
%    title(tm)
    

% Fourier Transforms 2222222222222222222222222222222222222222222222222
 figure(2)   
    pos = [0.33 0.1 0.3 0.5];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    hold on
    box on  
subplot(2,1,1)   % frequency
    xP = f; yP = PSD1;
    plot(xP,yP,'m','linewidth',2)
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
    legend(num2str(n1,'%2.0f'),num2str(n2,'%2.0f'),'orientation','horizontal','location','northeast')
    
subplot(2,1,2)    % wavelength
    xP = lambda; yP = PSD3;
    plot(xP,yP,'b','linewidth',2)
    xlabel('\lambda  ')
    ylabel('PSD [ a.u. ]')
  %  tm = sprintf('\\lambda_{Peak} = %3.2f     v = %3.2f   \n',max(wL_Peaks),v1);
    grid on
   % title(tm)
    set(gca,'fontsize',12)
    
% Numerical output ==================================================== 
% figure(4)   
%     pos = [0.65 0.1 0.3 0.6];
%     set(gcf,'Units','normalized');
%     set(gcf,'Position',pos);
%     hold on
%     box on    
%     
%    xlim([0 100])
%    ylim([0 160])
%    fs = 14;
%    
%   
%    Space = 160; space = 10;
%    txt = 'Model Parameters';
%    Htext = text(2,Space,txt);
%    set(Htext,'fontsize',fs,'color','k')
%   
%    Space = Space - space;
%       txt = ['   Atoms N = ' num2str(N,'%2.0f')];
%       Htext = text(2,Space,txt);
%       set(Htext,'fontsize',fs,'color','k')
%       
%    Space = Space - space; 
%       txt = ['  Time steps N_t = ' num2str(Nt,'%3.0f')];
%       Htext = text(2,Space,txt);
%       set(Htext,'fontsize',fs,'color','k')
%       
%    Space = Space - space; 
%       txt = ['   Simulation time t_{Max} = ' num2str(tMax,'%3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')
%       
%    Space = Space - space; 
%       txt = ['   d = ' num2str(d,'%3.2f') '    m = ' num2str(m,'%3.2f') '    k_S = ' num2str(kS,'%3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')
%       
%   Space = Space - space; 
%       txt = ['   \omega_{max} = ' num2str(w0,'%3.2f') '    f_{max} = ' num2str(f0,'%3.2f') '    v_{max} = ' num2str(v0,'%3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')
%       
%   Space = Space - space; 
%       txt = ['   Driving frequencies  f_1 = ' num2str(f1,'%3.2f') '    f_2 = ' num2str(f2,'%3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')    
%  
%   Space = Space - 2*space; 
%       txt = 'System response: Atom #11';
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')  
%       
%   Space = Space - space; 
%       txt = ['   f_{peak} = ' num2str(fPeaks1,'%3.2f    %3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')    
%      
%  Space = Space - space; 
%       txt = ['   \omega_{peak} = ' num2str(2.*pi.*fPeaks1,'%3.2f    %3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k') 
% 
%  if length(fPeaks1) ==  1   
%  Space = Space - space; 
%       txt = ['   Propagation speed v = ' num2str(v1,'%3.2f    %3.2f')  ];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')         
%  end
%  
%  Space = Space - space; 
%       txt = 'System response: Atom #41';
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')  
%       
%  Space = Space - space; 
%       txt = ['   f_{peak} = ' num2str(fPeaks2,'%3.2f    %3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')         
%       
%   Space = Space - space; 
%       txt = ['   \omega_{peak} = ' num2str(2.*pi.*fPeaks2,'%3.2f    %3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')      
%       
%  Space = Space - space; 
%       txt = 'Propagation parameters';
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')  
%       
%  Space = Space - space; 
%       txt = ['   \lambda_{peak} = ' num2str(wL_Peaks,'%3.2f    %3.2f')];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')     
%       
%  Space = Space - space; 
%       txt = ['   k_{peak} = ' num2str(2.*d./wL_Peaks,'%3.2f    %3.2f')  '  \pi/d'];
%       Htext = text(2,Space,txt); 
%       set(Htext,'fontsize',fs,'color','k')   
%       
%        
%       
%  axis off   
%    

    toc

 