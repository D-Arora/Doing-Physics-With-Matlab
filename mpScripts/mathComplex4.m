% mathComplex4.m

% Visualization of a Compound Complex Function / Fourier Component
% Phase Portrait at  a winding frequency fw

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% Matlab Version 2018a / 180915

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% John A Sims 
% email: john.sims@ufabc.edu.br
% Biomedical Engineering Department
% Federal University of ABC   Sao Bernardo Campus Brasil

clear
close all
clc

%  ANIMATED GIF =======================================================
   flagA = 1;       % save 1 / 0 not save
%  Enter file name
     ag_name = 'ag_A.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.00;  
%  Frame counter start
    nt = 1;

    
% INPUT SECTION =======================================================
% Winding frequency
fw = 1;

% Signal: amplitude R / frequency fs
%         number of periods numT / # elements nT
R = 1;
fs = 1; 
numT = 2;
nT = 251;


% CALCULATION SECTION  ================================================
  Ts = 1/fs;                    % signal period
  t = linspace(0,numT*Ts,nT);   % time elements

  
% Enter signal function ***********************************************
   % h = R.*sin(2*pi*fs*t) ;
   % h = R.*sin(2*pi*fs*t);
     h = 1.0.*sin(2*pi*fs*t) + 0.5.*sin(2*pi*2*fs*t) + 0.25.*sin(2*pi*3*fs*t);
% *********************************************************************

% Complex Exponential Function
   r = exp(-1j*2*pi*fw*t);
% Compound Complex Function   
   u = h .* r;
   uR = real(u);
   uI = imag(u);
   

% GRAPHICS SECTION: ANIMATED PHASE PORTRAIT ===========================   

%  Plot dimensions length
   L = max(abs(u));
   
figure(1)   
  pos = [0.05 0.05 0.30 0.35];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
    
  for cc = 1: nT
  % Unit Circle    
    tC = linspace(0,2*pi,200); xC = cos(2*pi*tC); yC = sin(2*pi*tC);
    plot(xC,yC)
    hold on
    xlim([-1.1*L,1.1*L])
    ylim([-1.1*L,1.1*L])
    
    grid on
    axis square
    xlabel('Re(u)','color','b') 
    ylabel('Im(u)','color','r')
    set(gca,'fontsize',14)
    box on
  
  % Phase portrait  
    xP = uR; yP = uI;
     plot(xP,yP,'k','linewidth',1)
    xP = [0 uR(cc)]; yP = [0 uI(cc)];
     plot(xP,yP,'k','linewidth',3)
    xP = [0 uR(cc)]; yP = [uI(cc) uI(cc)]; 
      plot(xP,yP,'b','linewidth',1.5)
    xP = [uR(cc) uR(cc)]; yP = [0 uI(cc)]; 
      plot(xP,yP,'r','linewidth',1.5)
    xP = uR(cc); yP = uI(cc); 
      Hplot = plot(xP,yP,'o');  
      set(Hplot,'markersize',8,'markerFaceColor','k','markerEdgeColor','k')
    xP = mean(uR); yP = mean(uI); 
      Hplot = plot(xP,yP,'o');  
      set(Hplot,'markersize',8,'markerFaceColor','m','markerEdgeColor','m')  
      
     
    tm1 = 'f_s = ';       tm2 = num2str(fs, '%2.1f'); tm3 = '  Hz';
    tm4 = '    f_w = ';  tm5 = num2str(fw,'%2.1f');
    tm6 = '  Hz';
    tm = [tm1 tm2 tm3 tm4 tm5 tm6];
    title(tm,'fontweight','normal')
      
    pause(0.01)
   
    hold off
   
% ANIMATION -------------------------------------------------------------
 if flagA == 1    
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
   hold off
  
  end
  
  
  
  
  
  
  
 