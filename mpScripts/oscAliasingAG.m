% oscAliasing.m
% Illustration of aliasing of the y-position signal using image that rotates at 1Hz
% Author: Sergio S Furuie
% Matlab 2018b  190417

% DOING PHYSICS ONLINE 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/oscAliasing.htm
% Ian Cooper  ian.cooper@sydney.edu.au
% School of Physics, University of Sydney



%%
% * Purpose *
%    To show visually and graphically the concept of aliasing in signal processing.
%    We are using a rotating spot with frequency 1 Hz and measuring the y-position y(t) of its center.
%    Therefore, y(t) is a sinusoidal signal with frequency 1Hz. 
%    The image is sampled and refreshed at different frame rates fs,showing visually the effects of sampling.
%    The signal y(t) is sampled with same sampling frequency fs.
%    We also indicate, in frequency axis, the 1Hz of signal and the Nyquist frequency, which is half of frame rate. 
%    When aliasing occurs, we also indicate the alias frequency, which is the actual signal frequency that will be discretized.
%    The alias frequency is a "reflection" of signal frequency in relation to Nyquist frequency.
%    Note what happens when alias frequency is less, equal and larger that 2xNyquist frequency.
%    In the latter case, the reflected frequency is negative and then is again reflected in relation to frequency 0.
%
% * Motivation *
%    For discrete signal and image data it is supposed that they were properly sampled.
%    Otherwise some degradations will be impossible to be corrected.
%    One of these degradations is aliasing, due to undersampling or lack of
%    analog filtering prior to sampling.
%
% * Fundaments * 
%     Let Y(f) be the Fourier transform of y(t).
%     It can be proved [1] that sampling y(t), with frequency fs, the resulting signal yp(t) has Fourier transform YP(f) that is a
%     repetition of a scaled Y(f) around multiples of fs: yP(f) = fs.sum{X(f-k.fs)};  k = -inf:inf.
%     Consider that Y(f) has significant components in the interval [-fm;fm]. Thus, the right part of Y(f-0.fs) will meet the left part of
%     Y(f-1.fs) at f = fs/2 when we increase fm. 
%     The frequency fs/2 is called Nyquist frequency of the discrete system. Any component beyond Nyquist frequency will overlap with 
%     other shifts of Y(f). In this case, even if we filter in the band [-fs/2;fs/2] the signal can not be recovered.
%     The Nyquist-Shannon sampling theorem says that if the highest meaningful component frequency of y(t) is fm, then the sampling
%     frequency should be at least 2.fm.
%
% * References *
%   1.Alan V. Oppenheim, Alan S. Willsky, S. Hamid, Signals and Systems, Pearson, ISBN-13: 978-0138147570

clear
clc
close all

% Animated gif
% Save animation as a gif file (0 no)  (1 yes) for flag2
    flagAG = 1;    
   

% Frequency of temporal signal (x-position) = rotation/s of spot  [Hz]
   f0 = 1;
% number of periods that is shown graphically   
   nCycles = 6;     

   msg = 'Click to proceed';

% Size of image
   Nx = 300;
   Ny = 300;

% Displaying rotated images accordingly to visualization rate
%   size of frame [left bottom width height]
   scrsz = get(groot,'ScreenSize');       
   pFull = 0.8;
   height = pFull*scrsz(4);
   b = scrsz(4)-height -scrsz(4)/10;
%   h = figure('Position',[10 b pFull*scrsz(3) height]);
   h = figure('Units','Normalized','Position',[0.05 0.2 0.4 0.6]);
   
%% [case: redundant] Nyquist freq >> Signal freq 
   frame1  = 0;
   frameRate = 8*f0;   % temporal sample frequency
   titulo = sprintf('Nyquist freq = %3.2f Hz  >>  Signal freq f_0 = %3.2f Hz   [case: redundant]',frameRate/2,f0);
   oscAliasingPlot(f0,frameRate,nCycles,Nx,Ny,titulo,h,flagAG,frame1);
%   uiwait(msgbox(msg,'Ok','modal'));

%% [case: fair]     Nyquist freq > Signal freq
   frame1 = frame1+1+fix(nCycles*frameRate/f0);
   frameRate = 3*f0;   % bom compromisso freq amostragem e quantidade de dados
   titulo = sprintf(' Nyquist freq = %3.2f Hz  >  Signal freq f_0 = %3.2f Hz   [case: fair]',frameRate/2,f0);
   oscAliasingPlot(f0,frameRate,nCycles,Nx,Ny,titulo,h,flagAG,frame1);
%   uiwait(msgbox(msg,'Ok','modal'));

%% [case: Nyquist criterion limit] Nyquist freq= Signal freq
   frame1 = frame1+1+fix(nCycles*frameRate/0);
   frameRate = 2*f0;   % [limite do criterio de Nyquist]
   titulo = sprintf(' Nyquist freq = %3.2f Hz  =  Signal freq f_0 = %3.2f Hz   [case: Nyquist criterion limit]',frameRate/2,f0);
   oscAliasingPlot(f0,frameRate,nCycles,Nx,Ny,titulo,h,flagAG,frame1);
%   uiwait(msgbox(msg,'Ok','modal'));

%% [case: aliasing < 2 Nyquist freq] Nyquist freq < Signal freq
   nCycles = 12;
   frame1 = frame1+1+fix(nCycles*frameRate/f0);   
   frameRate = 1.2*f0;   % aliasing
   titulo = sprintf('Nyquist freq = %3.2f Hz  <  Signal freq f_0 = %3.2f Hz   [case: aliasing < 2 Nyquist freq]',frameRate/2,f0);
   oscAliasingPlot(f0,frameRate,nCycles,Nx,Ny,titulo,h,flagAG,frame1);
%   uiwait(msgbox(msg,'Ok','modal'));

%% [case: aliasing = 2 Nyquist freq] Nyquist freq = (1/2) Signal freq
   frame1 = frame1+1+fix(nCycles*frameRate/f0);  
   frameRate = 1*f0;   % aliasing
   titulo = sprintf('Nyquist freq = %3.2f Hz  =  (1/2)Signal freq f_0 = %3.2f Hz  [case: aliasing = 2 Nyquist freq]',frameRate/2,f0);
   oscAliasingPlot(f0,frameRate,nCycles,Nx,Ny,titulo,h,flagAG,frame1);
%   uiwait(msgbox(msg,'Ok','modal'));

%% [case: aliasing > 2 Nyquist freq.] Nyquist freq. < (1/2)Signal freq.
   frame1 = frame1+1+fix(nCycles*frameRate/f0);
   frameRate = 0.8*f0;   % aliasing
   titulo = sprintf('Nyquist freq = %3.2f Hz  <  (1/2)Signal freq f_0 = %3.2f Hz   [case: aliasing > 2 Nyquist freq]',frameRate/2,f0);
   oscAliasingPlot(f0,frameRate,nCycles,Nx,Ny,titulo,h,flagAG,frame1);
%   msgbox('End of program','Ok');

%%
function oscAliasingPlot(freqSignal,frameRate,nCycles,Nx,Ny,titulo,h,flagAG,frame1)

% file name for animated gif   
    ag_name = 'ag_aliasing.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0.4; 

    
% Generate and show rotating spot, show x-position x(t) of spot center in relation to image center and frequencies
  clf(h)
  tFrame = 1/frameRate;
  nFrames = fix(nCycles*frameRate/freqSignal);  
  x_v = zeros(nFrames,1);
  dt_sig = 1/(50*freqSignal);             %sinal de referencia
  n_sig = nCycles/freqSignal/dt_sig;
  rot_init = pi;
  t_sig = (0:n_sig-1)*dt_sig;            %time values
  [~,~,~,dSpot_c] = oscAliasingSpot(Nx,Ny,rot_init);
  x_sig = dSpot_c*sin(rot_init-t_sig*freqSignal*2*pi);
  y_sig = dSpot_c*cos(rot_init-t_sig*freqSignal*2*pi);

  for n = 1:nFrames
    if(n==1),tic; end
    rot = rot_init-((n-1)*tFrame*freqSignal*2*pi);       %rotation (clockwise) of spot due to tFrame interval
    [circ_rot,xSpot,ySpot,dSpot_c] = oscAliasingSpot(Nx,Ny,rot);
    x_v(n) = xSpot;
    figure(h); 
    
    subplot(3,1,1)
    imagesc(1:Nx,1:Ny,circ_rot); axis image; xlabel('x'); ylabel('y');
    title(sprintf('Constant rotation f_0 = %3.2f Hz    Sampled freq f_S = %3.2f Hz',freqSignal,frameRate));
    t_v = (0:n-1)*tFrame; t = (n-1)*tFrame;
    set(gca,'fontsize',12)
    axis off
    
    subplot(3,1,2)
    hold on
    plot(t_sig,-y_sig,'k','linewidth',1);
    hPlot = plot(t_v,-x_v(1:n),'bo');
    set(hPlot,'markersize',6,'markerfacecolor','b')
    hPlot = plot(t,-xSpot,'mo');
    set(hPlot,'markersize',6,'markerfacecolor','m')
   % hPlot = plot(t_v,x_v(1:n),'bo',t,ySpot,'r*',t_sig,y_sig,'k:');
   % set(hPlot,'markersize',8,'markerfacecolor','b')
    grid on
    
    axis([0 (nCycles/freqSignal) -Nx/2 Nx/2]); 
    xlabel('time [s]'); ylabel('Y spot postion [a.u.]');  % 
    title(titulo);
    set(gca,'fontsize',12)
    box on
    
    % signal frequency and Nyquist frequency    
    fNyq = frameRate/2;
    
    subplot(3,1,3)
    line([freqSignal freqSignal],[0 1],'Color','b','LineWidth',3); line([fNyq fNyq], [0 0.7],'Color','r','LineWidth',3); 
    axis([0 5 0 1.2]);  xlabel('frequency [Hz]'); ylabel('a.u');
    set(gca,'fontsize',12)
    box on
    
    if(fNyq>=freqSignal)
        legend('Signal','Nyquist');
        title('Frequency domain: Signal freq and Nyquist freq');
    else      
        f_alias  = abs(2*fNyq - freqSignal);   % if may be larger than fNyq. Use other shifts of 2*fNyq
        f_aliasSign = sign(2*fNyq - freqSignal);
        while (f_alias > fNyq), f_alias = abs(2*fNyq-f_alias);f_aliasSign = sign(2*fNyq - f_alias); end
        line([f_alias f_alias],[0 1],'Color','g','LineWidth',3);
        legend('Signal','Nyquist','Alias');
        title(sprintf('Frequency domain: Signal freq, Nyquist freq and alias freq = %3.2f Hz',(-f_aliasSign)*f_alias));
    end 
    
    if flagAG > 0
         frame1 = frame1 + 1;
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
      % On the first loop, create the file. In subsequent loops, append.
         if frame1 == 1
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
         else
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
         end
    end
    
    drawnow;
    
     
    
 % Discounting processing time, if possible
    if(n==1), elapsedTime = toc; end
    tWait = max(tFrame-elapsedTime,0);
    pause(tWait);                   % discounting processing time
  end
  
       
end


% Building a circle with spot
function [circ,xSpot, ySpot,dSpot_c] = oscAliasingSpot(Nx,Ny,spotAngle)
% spotAngle in rad, counterclock.spotAngle=0 => x-axis 
% Using image notation: x-down, y-right

   R = fix(0.4*Nx);
   cx = fix(Nx/2)+1;
   cy = fix(Ny/2)+1;
   circ = zeros(Nx,Ny);
   spotR = fix(R/8);       %spot radius at (x,y)= (cx-R+spotR, 0)
   dSpot_c = R - spotR -2;  %distance of spot center to center of image
   xSpot = dSpot_c*cos(spotAngle);
   ySpot = dSpot_c*sin(spotAngle);
   cxSpot = cx + xSpot;
   cySpot = cy + ySpot;
   
for iy =1:Ny
    for ix =1:Nx
        d =sqrt((ix-cx)^2 +(iy-cy)^2);
        if (d>R), continue, end
        circ(ix,iy) = 1;
        dSpotCenter = sqrt((ix-cxSpot)^2 + (iy-cySpot)^2);
        if(dSpotCenter > spotR), continue, end
        circ(ix,iy) = circ(ix,iy) + 1;
    end
end
end

