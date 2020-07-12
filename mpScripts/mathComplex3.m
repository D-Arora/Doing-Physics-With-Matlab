% mathComplex3.m

% Visualization of a complex valued function
% The phase of the complex function is displayed using colors
% The calculation for the plots may take more than 100 seconds.

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


close all
clear
clc

N = 200;
x = linspace(-2,2,N);
y = x;


[xx, yy] = meshgrid(x,y);

zz = (xx+ 1j.*yy).^3 - 1 + 1j*3;
%zz = (xx+ 1j.*yy).^3 + 5;
zzR = real(zz); zzI = imag(zz);
zzA = abs(zz); zzP = angle(zz)./pi;


figure(1)
  pos = [0.02 0.5 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = xx; yP = yy; zP = zzR;
  surf(xP,yP,zP)
  view(31,18)
  shading interp
  colorbar
  set(gca,'xtick',-2:1:2)
  set(gca,'ytick',-2:1:2)
  grid on
  xlabel('x')
  ylabel('y')
  zlabel('z')
  title('Real(z)','fontweight','normal')
  set(gca,'fontsize',14)
  

  figure(2)
  pos = [0.33 0.5 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = xx; yP = yy; zP = zzI;
  surf(xP,yP,zP)
  view(19,28)
  shading interp
  colorbar
  set(gca,'xtick',-2:1:2)
  set(gca,'ytick',-2:1:2)
  grid on
  xlabel('x')
  ylabel('y')
  zlabel('z')
  title('Imag(z)','fontweight','normal')
  set(gca,'fontsize',14)
  
  
figure(3)
  pos = [0.65 0.5 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  sF = 1;   % scaling factor
  xP = xx; yP = yy; zP = zzA.^sF;
  surf(xP,yP,zP)
  view(19,28)
  shading interp
  colorbar
  set(gca,'xtick',-2:1:2)
  set(gca,'ytick',-2:1:2)
  grid on
  xlabel('x')
  ylabel('y')
  zlabel('z')
  title('scaled abs(z)','fontweight','normal')
  set(gca,'fontsize',14)

  
figure(4)
  pos = [0.02 0.05 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = xx; yP = yy; zP = zzP;
  surf(xP,yP,zP)
  view(22,56)
  shading interp
  HH = colorbar;
  set(HH,'Ticks',-1:0.25:1)
  set(HH,'Limits',[-1.1 1.1])
  set(gca,'xtick',-2:1:2)
  set(gca,'ytick',-2:1:2)
  grid on
  xlabel('x')
  ylabel('y')
  zlabel('z')
  title('phase(z)  [ rad/\pi ]','fontweight','normal')
  set(gca,'fontsize',14)
 
figure(5)
  pos = [0.33 0.05 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  s = 20;   % spacing of arrows
  P1 = xx(1:s:N,1:s:N); P2 = yy(1:s:N,1:s:N);
  P3 = zzR(1:s:N,1:s:N); P4 = zzI(1:s:N,1:s:N);
  HH = quiver(P1,P2,P3,P4);
  set(HH,'LineWidth',1,'MaxHeadSize',1,'AutoScaleFactor',1.5)
  set(gca,'xlim',[-2.5,2.5])
  set(gca,'ylim',[-2.5,2.5])
  axis square
  xlabel('x')
  ylabel('y')
  title('Vector Field','fontweight','normal')
  set(gca,'fontsize',14)
  set(gca,'fontsize',14) 
  
figure(6)
  pos = [0.65 0.05 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = xx; yP = yy; zP = zzA; v = 1:1:10;
 % contourf(xP,yP,zP,20,'ShowText','off');
  contourf(xP,yP,zP,v,'ShowText','on');
  shading interp
  colorbar
  hold on
  
  % plot unit circle with centre (0, 0)
  t = linspace(0,2*pi,200); R = sqrt(2);
  plot(R.*cos(t),R.*sin(t),'r','linewidth',2)
  
  set(gca,'xtick',-2:1:2)
  set(gca,'ytick',-2:1:2)
  grid on
  xlabel('x')
  ylabel('y')
  zlabel('z')
  title('abs(z)','fontweight','normal')
  axis square
  set(gca,'fontsize',14) 