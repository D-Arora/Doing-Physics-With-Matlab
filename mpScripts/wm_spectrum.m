% wm_spectrum.m

% Color plot of the visible spectrum
%   for wavelengths from 370 nm to 770 nm
% Calls the function ColorCode.m

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180712

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm

clear 
close all
clc

figure(1)
    pos = [0.1 0.1 0.3 0.2];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    set(gcf,'color','w');
    N = 512;
    xP = linspace(380,780,N);
    yP = ones(1,length(xP));
    
    hold on
    
    thisColorMap = hsv(512);
    
for cx = 1:N-1
    wL = xP(cx)*1e-9;
    thisColor = ColorCode(wL);
    h_area = area(xP(cx:cx+1),yP(cx:cx+1));
    set(h_area,'FaceColor',thisColor);
    set(h_area,'EdgeColor',thisColor); 
    set(gca,'xLim',[380 770]);
end
       
    xlabel('wavelength  \lambda   [ nm ] ','fontsize',14);
    set(gca,'fontsize',14);
    set(gca,'yTick',[]);
       
 
 