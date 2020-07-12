% em_swe_01.m

% 14 nov 15
% Solving the scalar wave equation using the FDTD method
% Ian Cooper
% School of Physics, University of Sydney
% documentation: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts

close all
clear all
clc
tic
% ========================================================================
% INPUTS
% ========================================================================

% Number of spatial grid points
  Nx = 100;
% Number of time steps
  Nt = 25;
% length of simulation region
   L = 100;
% wave or propagation speed
   v = 10;
% Courant number
   S = 1;

% ========================================================================
% SETUP
% ========================================================================

% Spatial grid spacing / time step / spatial grid
   hx = L / Nx;
   ht = S * hx / v;
   S2 = S^2;
   x = linspace(0, L , Nx);
   t = linspace(0,ht*Nt,Nt);
% Initialise wavefunction   
   u = zeros(Nx,Nt);

% Initialize wavefunction - pulse type (comment / uncomment flag)
    flag = 1;    % Gaussian pulse
   %  flag = 2;    % Rectangular pulse

if flag == 1
% Gaussian pulse inputs / 1st time step / 2nd time step
   A = 1;        % pulse amplitude
   s = 2;        % pulse width
   np = 80;      % pulse starting point given by index np
for nx = 2 : Nx-1   % time step nt = 1
   %  u(nx,1) = A .* exp(-0.5 .* ((x(nx)- x(np))./(s*hx)).^2);
    u(nx,1) = A .* exp(-0.5 .* ((nx-np)./s).^2);
end

for nx = 2 : Nx - 1   % time step 2
   %  u(nx,2) = A .* exp(-0.5 .* ((x(nx)-x(np) - v*ht)./(s*hx)).^2);
      u(nx,2) = A .* exp(-0.5 .* ((nx-np-1)./s).^2);
     %u(nx,2) = u(nx,1) - 0.5*S2*(u(nx+1,1) - 2*u(nx,1) + u(nx-1,1));
     %u(nx,2) =2*u(nx,1)-u(nx,1) + S2*(u(nx+1,1) - 2*u(nx,1)+u(nx-1,1)) ;  
 end
  %u(:,2) = u(:,1);
end

if flag == 2
% Rectangular pulse
   u(10:30,1) = 0.5;
   u(11:31,2) = 0.5;
end
%
%
% ========================================================================
% Solving the scalar wave equation FDTD   
% ========================================================================

for nt = 2 : Nt-1
    for nx = 2 : Nx-1
        u(nx,nt+1) = 2*u(nx,nt)- u(nx,nt-1)+ S2*(u(nx+1,nt) - 2*u(nx,nt) + u(nx-1,nt));     
    end
  % Boundary Conditions    comment / uncomment
      u(1, nt+1)  = 0;               % Fixed end
     u(Nx, nt+1) = 0;               % Fixed end
    %  u(Nx,nt+1) = u(Nx-1,nt+1);      % Free end
     % u(1,nt+1) = u(2,nt+1);          % Free end
    %  u(Nx,nt+1) = u(Nx-1,nt);      % ABC
end

% ========================================================================
% GRAPHICS
% ========================================================================
 
   figure(2)
      set(gcf,'units','normalized','position',[0.5,0.1,0.4,0.7]); 
      grid on
      col = 'b'; LW = 2; fs = 14;
      tm1 = 't  =  ';
      tm3 = '  s ';
      subplot(5,1,1)
         xP = x; yP = u(:,1);
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1])
         tm2 = num2str(t(1),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);
      subplot(5,1,2)
         xP = x; yP = u(:,round(Nt/4));
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1])
         tm2 = num2str(t(round(1*Nt/4)),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);
     subplot(5,1,3)
         xP = x; yP = u(:,round(2*Nt/4));
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1]);
         tm2 = num2str(t(round(2*Nt/4)),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);
     subplot(5,1,4)
         xP = x; yP = u(:,round(3*Nt/4));
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1]);
         tm2 = num2str(t(round(3*Nt/4)),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);
     subplot(5,1,5)
         xP = x; yP = u(:,round(Nt));
         plot(xP, yP,col,'linewidth',LW);
         axis([0 L -1.1 1.1]);
         tm2 = num2str(t(Nt),'%2.1f\n');
         tm = [tm1 tm2 tm3];
         title(tm);    
 
 % ========================================================================        
 %   ANIMATION: SAVE FILE AS animated gig
 % ========================================================================
 
   f_gif = 1;                 % flag: 0 animated gif NOT saved or 1 saved
   ag_name = 'ag_em_001.gif'; % file name for animated gif
   delay = 0.1;  % A scalar value 0 to 655 (inclusive)- specifies
                 % the delay in seconds before displaying the next image
 
% Setup for saving images (im)

 figure(1)
    set(gcf,'units','normalized','position',[0.05,0.1,0.4,0.3]);
    set(gcf,'color','w');
    grid on
    col = 'b'; LW = 2; fs = 14;
    c = 1;
    tm1 = 't  =  ';
    tm3 = '  s ';
  
    xP = x; yP = u(:,1); col = [0.5 0.5 0.5]; LW = 1;
    plot(xP, yP,'color',col,'linewidth',LW);
    axis([0 L -1.1 1.1]);
    
   % Ct = Nt/10;
    f = getframe(gcf);
    [im,map] = rgb2ind(f.cdata,256,'nodither');  %RGB to indexed images
        im(1,1,1,Nt) = 0;  
    
    
   LW = 2; col = [0.2 0.2 1];
   for nt = 1: Nt
      xP = x; yP = u(:,nt);
      plot(xP,yP,'color',col,'linewidth',LW);
      axis([0 L -1.1 1.1])
      grid on
      set(gca,'fontsize',fs);
      xlabel('x'); ylabel('u');
      tm2 = num2str(t(nt),3);
      tm = [tm1 tm2 tm3];
      title(tm);
      pause(0.1);
      drawnow;
     f = getframe(gcf);
     im(:,:,1,c) = rgb2ind(f.cdata,map,'nodither');
     c = c+1;
   end
   
     
         
if f_gif > 0

% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously

imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);

end     % end if 
% ========================================================================
toc