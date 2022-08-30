% cemEDR02.m

% ELECTRIC DIPOLE RADIATION
%   GENERATION OF ELECTROMAGNETIC WAVES
%   Dipole in Z direction

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 220815 / Matlab version R2021b

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Notes
%    https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cemEDR.htm
% SCRIPTS
%    https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb



close all
clc
clear

% Animated gif  =======================================================
% Show and Save animation as a gif file (0 no)  (1 yes) for flag2
    flagAG = 0;    
% file name for animated gif   
    ag_name = 'ag_cemEDR02.gif'; 
% Delay in seconds before displaying the next image  
    delay = 0;
    frame1 = 0;

constantsEM
% speed of light c
% Coulomb constant ke

% time 
  Nt = 99; tMin = 0; tMax = 20e-9;
  t = linspace(tMin,tMax, Nt);

% Polar angle  
  theta = 1*pi/2;

% wavelength wL / propagation constant k / freq f / ang freq  w
  wL =  1;
  k = 2*pi/wL;
  f = c/wL;
  T = 1/f;
  w = 2*pi*f;
  wt = w.*t;

% Radial displacement
  N = 1299; rMin = 5; rMax = 25; %25;
  r = linspace(rMin,rMax,N);
  kr = k.*r;

% Compute Electric field: polar component 
%         Magnetic field: azimuthal component
%         Poynting vector
  EP = zeros(N,Nt); BA = zeros(N,Nt); ER = zeros(N,Nt);
  for nt = 1 : Nt
      
       ER(:,nt) = 10 .* (1 + 1i./kr).* exp(1i.*(kr - wt(nt)))./r.^2;
       EP(:,nt) = 1.*(sin(theta)./r) .* ( -sin(kr - wt(nt))./kr.^2 + cos(kr - wt(nt))./kr + sin(kr - wt(nt)) );
       % EP(:,nt) = 1.* sin(theta).*exp(1i.*(kr - wt(nt)))./(r); 
     %  BA(:,nt) =  10.*(sin(theta)./r) .* ((1./kr) .* cos(kr - wt(nt)) - sin(kr - wt(nt)));
        BA(:,nt) = (sin(theta)./r).*( cos(kr - wt(nt))./kr - sin(kr - wt(nt)));
      % EP(:,nt) = 10.*(sin(theta)./r) .* ( (1./kr) .* cos(kr - wt(nt)) + (1 - 1./(kr).^2) .* sin(kr - wt(nt))) ;
      % EP(:,nt) = 10.*(sin(theta)./r.^2) .* (sin(kr - wt(nt)) - kr.*cos(kr - wt(nt)));
  end
    % EP = real(EP); ER = real(ER); BA = real(BA);
  Savg = sin(theta./r).^2;


% GRAPHICS  =======================================================
figure(1)
    set(gcf,'Units','normalized');
    set(gcf,'Position',[0.1 0.1 0.30,0.45])

 for nt = 1:Nt
subplot(2,1,2)
   xP = r; yP = BA(:,nt);
      plot(xP,yP,'k','linewidth',2)
   hold on
   xP = r(200); yP = BA(200,nt);
      Hplot = plot(xP,yP,'ok','linewidth',1);
   set(Hplot,'markerfacecolor','k')
  xP = r; yP = 5.*max(max(BA))./(r);
    plot(xP,yP,'r','linewidth',1)

 % ylim([-2 2])
   xlim([rMin rMax])
  grid on; box on
  set(gca,'FontSize',14)
  xlabel('r  [m]'); ylabel('B_\phi field  [a.u.]')
  pause(0.01)
  hold off

subplot(2,1,1)
  xP = r; yP = EP(:,nt);
    plot(xP,yP,'b','linewidth',2)
  hold on
  xP = r; yP = 5.*max(max(EP))./r;
    plot(xP,yP,'r','linewidth',1)
  xP = r(200); yP = EP(200,nt);
    Hplot = plot(xP,yP,'ob','linewidth',1);
    set(Hplot,'markerfacecolor','b')
 % ylim([-2.5 2.5])
  xlim([rMin rMax])
  grid on; box on
  set(gca,'FontSize',14)
%  legend('E_\theta',' ','E_r','|E|','B_\phi','Orientation','horizontal','Location','bestoutside')
  xlabel('r  [m]'); ylabel('E_\theta field  [a.u.]')
  txt = sprintf('\\theta = %2.3f  rad/\\pi  \n',theta/pi);
  title(txt)
  pause(0.01)

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
  hold off


end
%   =====================================================================
figure(2)
   set(gcf,'Units','normalized');
   set(gcf,'Position',[0.42 0.1 0.25,0.25])
   xP = t.*1e9; yP = EP(100,:);
   plot(xP,yP,'b','linewidth',2)
   grid on; box on
   xlabel('t  [ns]'); ylabel('E_\theta  [a.u.]')
   set(gca,'FontSize',14)

   %%
figure(3)
 clf(3)
  set(gcf,'Units','normalized');
  set(gcf,'Position',[0.68 0.1 0.30,0.25])
  thetaA = linspace(-pi,pi,199);
  R = [0.5 0.6 0.7 0.8 0.9 1];
   
  for n = 1:length(R)
  rho = (sin(thetaA)./R(n)).^2;
    hh = polarplot(thetaA,rho,'linewidth',2);
    hold on
  end
  title('S_{avg}')
  rticks([ ])
  rlim([0 4])
  pax = gca;
  pax.ThetaDir = 'clockwise';
  pax.ThetaZeroLocation = 'top';
% pax = gca;
%  pax.ThetaAxisUnits = 'radians';
 legend('0.5','0.6','0.7','0.8','0.9','1','location','best')
  grid on
  set(gca,'FontSize',14)

%%
% close
%   figure(4)
% 
%   set(gcf,'Units','normalized');
%   set(gcf,'Position',[0.68 0.4 0.30,0.4])
% R=5; % outer radius of torus
% r=4.8; % inner tube radius
% th=linspace(0,2*pi,99); % e.g. 36 partitions along perimeter of the tube 
% phi=linspace(0,2*pi,99); % e.g. 18 partitions along azimuth of torus
% % we convert our vectors phi and th to [n x n] matrices with meshgrid command:
% [Phi,Th]=meshgrid(phi,th); 
% % now we generate n x n matrices for x,y,z according to eqn of torus
% x=(R+r.*cos(Th)).*cos(Phi);
% y=(R+r.*cos(Th)).*sin(Phi);
% z=sin(Th).^2;
% surf(x,y,z); % plot surface
% hold on
% surf(x,y,-z); % plot surface
% daspect([1 1 0.2]) % preserves the shape of torus
% colormap('winter') % change color appearance 
% shading interp
% title('Torus')
% xlabel('X');ylabel('Y');zlabel('Z');
% view(-15,42)
% box on

%%
%clear; clc; close all
  figure(4)

  set(gcf,'Units','normalized');
  set(gcf,'Position',[0.68 0.4 0.30,0.4])
%r = @(u,v) 2 + sin(7.*u + 5.*v);
 r = @(u,v) sin(v).^2;
funx = @(u,v) r(u,v).*cos(u).*sin(v);
funy = @(u,v) r(u,v).*sin(u).*sin(v);
funz = @(u,v) r(u,v).*cos(v);
fsurf(funx,funy,funz,[0 2*pi 0 pi]) 
%camlightr = @(u,v) 2 + sin(7.*u + 5.*v);
funx = @(u,v) r(u,v).*cos(u).*sin(v);
funy = @(u,v) r(u,v).*sin(u).*sin(v);
funz = @(u,v) r(u,v).*cos(v);
fsurf(funx,funy,funz,[0 2*pi 0 pi]) 
view(-40,50)
axis off
%camlight
