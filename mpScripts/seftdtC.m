% sefdtdC.m


clear
close all
clc

tic
Nt = 10000;
L = 4e-9;
Nx = 1001;

C1 = 0.1;

% Initialise the pulse
   nx0 = round(Nx/4);     % pulse centre   round(Nx/4)
   %wL = L/40;             % pulse wavelength   L/40
   s = L/25;              % pulse width        L/25
   wL = 1.6e-10;

x = linspace(0,L,Nx); dx = x(2) - x(1);

% for nx = 1 : Nx
%    yR(nx) = exp(-0.5*((x(nx)-x(nx0))/s)^2)*cos(2*pi*(x(nx)-x(nx0))/wL);
%    yI(nx) = exp(-0.5*((x(nx)-x(nx0))/s)^2)*sin(2*pi*(x(nx)-x(nx0))/wL);
% end

   yR = exp(-0.5.*((x-x(nx0))./s).^2).*cos((2*pi).*(x-x(nx0))./wL);
   yI = exp(-0.5.*((x-x(nx0))./s).^2).*sin((2*pi).*(x-x(nx0))./wL);
   
   psiR(1,:) = yR;
   psiI(1,:) = yI;
   
  C2 = 0;
  U = zeros(1,Nx);
  psiR = zeros(Nt,Nx);
  psiI = zeros(Nt,Nx);

for nt = 1 : Nt-1
   for nx = 2 : Nx - 1
      yR(nx) = yR(nx) - C1*(yI(nx+1)-2*yI(nx)+yI(nx-1)) + C2*U(nx)*yI(nx);
    end
   psiR(nt+1,:) = yR;
   
   for nx = 2 : Nx-1
      yI(nx) = yI(nx) + C1*(yR(nx+1)-2*yR(nx)+yR(nx-1)) - C2*U(nx)*yR(nx);
   end
    psiI(nt+1,:) = yI;
end 



figure(1)
  pos = [0.05 0.05 0.3 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  FS = 14;
  xP = x;
  
  subplot(2,1,1)
  yP = yR;
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = psiR(Nt,:);
  plot(xP,yP,'r','linewidth',2)
  grid on
  box on
  set(gca,'fontsize',FS)
  
  
%    tm1 = 't_{steps} = ';
%    tm2 = num2str(nt,'%5.0f\n');
%    tm3 = '  <E> = ';
%    tm4 = num2str(real(Eavg),'%3.0f\n');
%    tm5 = '  <K> = ';
%    tm6 = num2str(real(Kavg),'%3.0f\n');
%    tm7 = '  <U> = ';
%    tm8 = num2str(real(Uavg),'%3.0f\n');
%    tm9 = '   eV';
%    tm = [tm1 tm2 tm3 tm4 tm5 tm6 tm7 tm8 tm9];
%    h_title = title(tm,'fontsize',fs);
%   
  subplot(2,1,2)
  yP = yI;
  plot(xP,yP,'b','linewidth',2)
  hold on
  yP = psiI(Nt,:);
  plot(xP,yP,'r','linewidth',2)
  grid on
  box on
  set(gca,'fontsize',FS)
  
 % ylim([-1e5 1e5])
  toc
  
 