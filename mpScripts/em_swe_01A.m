% em_swe_01.m

% 14 nov 15
% Solving the scalar wave equation using the FDTD method
% Ian Cooper
% School of Physics, University of Sydney
% documentation: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts

close all
clear 
clc
tic

% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flagS = 1 (save)
   flagS = 1;
%  Enter file name
     ag_name = 'ag_A.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.1;  
%  Frame counter start
%    nt = 1;


% ========================================================================
% INPUTS
% ========================================================================

% Number of spatial grid points
  Nx = 301;
% Number of time steps
  Nt = 2201;
% length of simulation region
   L = 10;
% wave or propagation speed
   v = 20;
% Courant number
   S = 1;
% Thee signal detection location SL / Clamped string location CL  (indices) 
  SL = round([Nx/2, Nx/3,4*Nx/5]);
 % CL = round([Nx/3 2*Nx/3]);
   CL = round(Nx/2);
 %  CL = 1;
 % Driving frequency  >>>>
   f0  = 1;
% ========================================================================
% SETUP
% ========================================================================

% Spatial grid spacing / time step / spatial grid
   hx = L / Nx;
   ht = S * hx / v;
   S2 = S^2;
   x = linspace(0, L , Nx);
   t = linspace(0,ht*Nt,Nt);

% Random initial conditions  (comment / uncommment selection)
   u = zeros(Nx,Nt);
   u(2:Nx-1,1) = 1 .* (2.*rand(Nx-2,1) - 1);
%  u(2:Nx-1,1) = 0.25.*rand(Nx-2,1);
%   u(2:51,1) = (1/5).*x(2:51) - 5.1/5;
%    u(52:Nx,1) = -(1/5).*x(52:Nx) + 5.5/5;
%   u(38:40,1) = 0.5;

% ========================================================================
% Solving the scalar wave equation FDTD   
% ========================================================================

for nt = 2 : Nt-1
    for nx = 2 : Nx-1
    % sinusoidal excited at x = 0;
      u(1,nt+1) = 0.2*sin(2*pi*f0*t(nt));
        u(nx,nt+1) = 2*u(nx,nt)- u(nx,nt-1)+ S2*(u(nx+1,nt) - 2*u(nx,nt) + u(nx-1,nt));
% Clamped location   (comment / uncomment     
%        u(CL,nt+1) = 0;

    end
% Boundary Conditions    comment / uncomment
      u(1,  nt)  = 0;               % Fixed end
      u(Nx, nt) = 0;               % Fixed end
end

     
% ========================================================================
% FOURIER TRANSFORMS
% ========================================================================
% Frequency range
  fMax = 10;
  f = linspace(0,fMax,Nt);
% Initialize F.T.  
%  HB = zeros(1,Nt); HR = zeros(1,Nt); HM = zeros(1,Nt);

  HB = FT(u(SL(1),:),f,t,Nt);
  psdB = HB.*conj(HB);

  HR = FT(u(SL(2),:),f,t,Nt);
  psdR = HR.*conj(HR);
 
  HM = FT(u(SL(3),:),f,t,Nt);
  psdM = HM.*conj(HM);


% ========================================================================
% GRAPHICS
% ========================================================================
 
figure(2)
      set(gcf,'units','normalized','position',[0.5,0.1,0.25,0.55]); 
     
subplot(3,1,1)
      xP = t; yP = u(SL(1),:);
      plot(xP,yP,'b','linewidth',2)
      hold on
      yP = u(SL(2),:);
      plot(xP,yP,'r','linewidth',2)
      yP = u(SL(3),:);
      plot(xP,yP,'m','linewidth',2)
      xlabel('t  [s]'); ylabel('u  [m]')
      set(gca,'FontSize',14)
      grid on; box on

subplot(3,1,2)
      xP = f; yP = psdB;
      plot(xP,yP,'b','linewidth',2)
      hold on
      yP = psdR;
      plot(xP,yP,'r','linewidth',2)
       yP = psdM;
      plot(xP,yP,'m','linewidth',2)
      xlabel('f  [Hz]'); ylabel('psd')
      set(gca,'FontSize',14)
      grid on; box on

subplot(3,1,3)
      xP = f; yP = log(psdB);
      plot(xP,yP,'b','linewidth',2)
      hold on
      yP = log(psdR);
      plot(xP,yP,'r','linewidth',2)
      yP = log(psdM);
      plot(xP,yP,'m','linewidth',2)
      ylim([-20 10]); 
      grid on; box on
      xlabel('f  [Hz]'); ylabel('log(psd)')
      set(gca,'FontSize',14)
 

 figure(1)  % 11111111111111111111111111111111111111111111111111111111
    set(gcf,'units','normalized','position',[0.05,0.1,0.3,0.25]);
    set(gcf,'color','w');
    grid on
    c = 1;
    tm1 = 't  =  ';
    tm3 = '  s ';
  
    xP = x; yP = u(:,1); col = [0.5 0.5 0.5]; LW = 1;
    plot(xP, yP,'color',col,'linewidth',LW);
    axis([0 L -10.1 10.1]);
     

   for nt = 1:5: Nt
      xP = x; yP = u(:,nt);
      plot(xP,yP,'k','linewidth',2);
      hold on
      xP = x(SL(1)); yP = u(SL(1),nt);
      plot(xP,yP,'bo','MarkerFaceColor','b','MarkerSize',10)
      xP = x(SL(2)); yP = u(SL(2),nt);
      plot(xP,yP,'ro','MarkerFaceColor','r','MarkerSize',10)
      xP = x(SL(3)); yP = u(SL(3),nt);
      plot(xP,yP,'mo','MarkerFaceColor','m','MarkerSize',10)
      %axis([0 L -20.1 20.1])
      ylim([-1.1*max(max(u)); 1.1*max(max(u))])
      grid on
      set(gca,'fontsize',14);
      xlabel('x'); ylabel('u');
      txt = sprintf('t = %2.2f s   f0 = %2.2f Hz  \n   ',t(nt),f0);
      tm2 = num2str(t(nt),3);
      tm = [tm1 tm2 tm3];
      title(txt);

  if flagS == 1
       frame = getframe(1);
       im = frame2im(frame);
       [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
       if nt == 1
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
       else
         imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
       end
    %   nt = nt+1;
   end

      pause(0.001);
      hold off
   end
   
             

% ========================================================================

toc

function  H = FT(u1,f,t,Nt)
  H = zeros(Nt,1);
 % Fourier Transform H(f)
  for c = 1:Nt
   g = u1 .* exp(1i*2*pi*f(c)*t);
   H(c) = simpson1d(g,min(t),max(t));
 end


end
