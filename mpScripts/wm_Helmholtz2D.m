% wm_Helmholtz2D.m

% Helmholtz Equation: Solving the eigenvalue problem
%                     for transverse standing waves on a square membrane

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/wm_Helmholtz2D.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  181104

close all
clear 
clc
tic

% INPUTS  =============================================================

% Mode number for graphical outputs  
mX = 2; mY = 3;
% Number of grid points
N = 59;
% Number of time steps
nT = 100;
% Length of rod [m]
LX = 1; LY = 2;
% speed of transverse wave along rod [m/s]
v = 300;


% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flagS = 1 (save)
   flagS = 1;
%  Enter file name
     ag_name = 'ag_A.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.1;  
%  Frame counter start
    nt = 1;


% CALCULATIONS  ======================================================
% Spatial domian  [m]
  x = (0:N+1).*(LX/(N+1));
  y = (0:N+1).*(LY/(N+1));
  
% Eigenvalue Matrix A: eigenfunctions (eignFN) / eigenvalues (eignV) 
  off = ones(N-1,1);
  A = 2*eye(N) - diag(off,1) - diag(off,-1);

  [eignFN, eignV] = eig(A);
  
% Spatial Wavefunction  US
  USX = zeros(N+2,1);
  USX(2:N+1) = eignFN(:,mX);
  USX = USX ./max(USX);

  USY = zeros(N+2,1);
  USY(2:N+1) = eignFN(:,mY);
  USY = USY ./max(USY);
  
  USXX = meshgrid(USX);
  USYY = meshgrid(USY)';
  
  USS = USXX.*USYY;
  
% Time dependent wavefunction UT
  % propagation constant  [1/m]
    kX = sqrt(eignV(mX,mX)) .* (N+1)/LX;
    kY = sqrt(eignV(mY,mY)) .* (N+1)/LY;
    k = sqrt(kX^2 + kY^2);
  % angular frequency  [rad/s]
    w = v*k;
  % period  [s]
    T = 2*pi/w;
  % frequency  [Hz]
    f = 1/T;
  % time  [s]
    t = linspace(0,1*T,nT);
  % wavelength [m]
    lambda = 2*pi/k;
 % time dependent wavefunction   
    UT = cos(w.*t);   
  
  
% GRAPHICS  ==========================================================
 figure(1)
  FS = 12;
  pos = [0.05 0.05 0.40 0.4];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
 
for cc = 1:nT
  subplot('Position', [0.05 0.4 0.9 0.6])
    mesh(x,x,USS.*UT(cc))
    shading interp 
    colormap(summer)
    view(-50,40)
    zlim([-1.1 1.1])
    axis off
  
  subplot('Position', [0.05 0.02 0.9 0.3])
    zP = USS.*UT(cc); zP(1,1) = -1; zP(N,N);
    contourf(x,y,zP)
    shading interp 
    colormap(summer)
    view(-50,90)
    zlim([-1.1 1.1])
    axis off
    h = colorbar('limits',[-1 1]);
    tm1 = 'Mode   m_X =  ';
    tm2 = num2str(mX,'%3.0f     m_Y =     \n'); 
    tm3 = '  ';
    tm4 = num2str(mY,'%5.0f \n');
    tm5 = '     f =  ';
    tm6 = num2str(f,'%5.0f  Hz \n');
    tm = [tm1 tm2 tm3 tm4 tm5 tm6];
    title(tm,'fontsize',FS)
    
          
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
       nt = nt+1;
   end
    
   pause(0.1)
end

figure(2)
  FS = 12;
  pos = [0.6 0.05 0.3 0.2];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = x;  yP = USX;
    plot(xP,yP,'b','linewidth',2);
  hold on
  xP = y;  yP = USY;
    plot(xP,yP,'r','linewidth',2);
  grid on
  set(gca,'fontsize',FS)
  xlabel('position along X & Y axes')
  ylabel('wavefunctions')
  legend('X','Y','location','southeast')
  
    tm1 = 'Mode   m_X =  ';
    tm2 = num2str(mX,'%3.0f     m_Y =     \n'); 
    tm3 = '  ';
    tm4 = num2str(mY,'%5.0f \n');
    tm5 = '     f =  ';
    tm6 = num2str(f,'%5.0f  Hz \n');
    tm = [tm1 tm2 tm3 tm4 tm5 tm6];
    title(tm,'fontsize',FS)
  
toc
