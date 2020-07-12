% cemPol1.m

% 180118
% Polarization of Light 
% Plot resultant electric field from two sinusoidally time varying
%   orthogonal components.  
% We can observe linear, eliptical and circular polarization. 
% We can also observe the right or left handedness of the polarization for
%   the ellipcal and circular cases as viewed by an observer 
%   looking towards the source of the radiation.

% John A Sims 
% email: john.sims@ufabc.edu.br
% Universidade Federal do ABC 
% http://www.ufabc.edu.br/
% Biomedical Engineering Department,
% Centro de Engenharia, Modelagem e Ciencias Sociais Aplicadas.
% (Centre for Engineering, modelling and applied social sciences)
% Rua Arcturus (Jd Antares)
% Anchieta
% 09606070 - Sao Bernardo do Campo, SP - Brasil

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm


close all
clear all
clc

% =======================================================================
%  INPUTS  default value [  ] 
% =======================================================================
%  Electric field amplitudes [10] [10]
    E0x =  10; E0y = 10; 
% Phase of Ey in radians  [0.5*pi]
%  Try phi [rad]: 0.25*pi 0.5*pi 0.75*pi pi 0 -0.25*pi -0.5*pi -0.75*pi -pi 
    phi = 0.5*pi;
% Grid points
   N = 24;

% =======================================================================   
% Setup for saving images animated gif file and AVI files
% =======================================================================
% Save animation flagA = 0 or 1 or 2  [default flagA = 0]
   flagA = 0;
   % flagA = 0   Animation NOT saved
   % flagA = 1   Animation saved as an animated gif file
   % flagA = 2   Animation saved as an AVI file

switch flagA
case 1
%  ANIMATED GIF:
     ag_name = 'ag_Pol01.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.20;  
%  Frame counter start
    nt = 1;

case 2   
%  AVI: open writer object for avi.
     aviName = 'avi_Pol01.avi';
     aviObj = VideoWriter(aviName);
     aviObj.FrameRate = 10;  % from 30
     aviObj.Quality = 40;    % from 75
     open(aviObj);
end   %  switch   
   
   
% CALCULATIONS z = 0 =====================================================
   wt = linspace(0,2*pi,N);
% Electric fields
   Ex = E0x .* exp(1j*wt);
   Ey = E0y .* exp(1j*(wt + phi));
   E = real(Ex) + 1j.*real(Ey);

% GRAPHICS  =============================================================
 
 figure(1)   
    pos = [0.07 0.05 0.25 0.35];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    set(gcf,'color','w');
    fss = 12;
    
 for cc = 1 : N
    xP = real(E); yP = imag(E); 
    plot(xP, yP,'b','linewidth',0.5); 
    set(gca,'xLim',[-12 12]);
    set(gca,'yLim',[-12 12]);
    grid on
    xlabel('E_x')
    ylabel('E_y');
    tm1 = '\phi = ( ';
    tm2 = num2str(phi / pi,  '%3.2f \n');
    tm3 = ' ) \pi  rad ';
    tm = [tm1 tm2 tm3];
    title(tm,'fontweight','normal');
    axis square
    set(gca,'fontsize',fss)
    
    hold on
    for c = 1: N
      xP = [0 real(E(c))]; yP = [0 imag(E(c))]; 
      plot(xP, yP,'col',[0.5 0.5 0.5],'linewidth',0.2)
    end
   
   xP = real(E(cc)); yP = imag(E(cc));
   arrow([0,0], [xP,yP],'Width',2, 'EdgeColor',[1 0 0],'FaceColor',[1 0 0 ]);
   pause(0.2)
   hold off

  % Adding the frame to the video object using the 'writeVideo' method
   if flagA == 2
     frame = getframe(1);
     writeVideo(aviObj,frame);
   end 
  
  % Adding the frame to the animated gif 
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
 end
    
 if flagA == 2
  close(aviObj);
end
