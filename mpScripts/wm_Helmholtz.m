% wm_Helmholtz.m

% Helmholtz Equation: Solving the eigenvalue problem
%                     for transverse standing waves on a rod

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/wm_Helmholtz.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  181104

close all
clear 
clc
tic

% INPUTS  =============================================================
% Boundary conditions:
%     fixed fixed BC = 1 / fixed free BC = 2 / free free BC = 3 
BC = 2; 
% Mode number for graphical outputs  BC = 3 --> m > 1
m = 1;
% Number of grid points
N = 599;
% Number of time steps
nT = 100;
% Length of rod [m]
L = 1;
% speed of transverse wave along rod [m/s]
v = 300;


% =====================================================================
%  ANIMATED GIF:   flagS = 0 (not saved)  / flagS = 1 (save)
   flagS = 0;
%  Enter file name
     ag_name = 'ag_A.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.1;  
%  Frame counter start
    nt = 1;


% CALCULATIONS  ======================================================
% Spatial domian  [m]
  x = (0:N+1).*(L/(N+1));

% Eigenvalue Matrix A: eigenfunctions (eignFN) / eigenvalues (eignV) 
  off = ones(N-1,1);
  A = 2*eye(N) - diag(off,1) - diag(off,-1);

  if BC == 2; A(N,N) = 1; end               % fixed free
  if BC == 3; A(1,1) = 1; A(N,N) = 1; end   % free free

  [eignFN, eignV] = eig(A);
  
% Spatial Wavefunction  US
  US = zeros(N+2,1);
  US(2:N+1) = eignFN(:,m);
  US = US ./max(US);

  if BC == 2; US(N+2) = US(N+1); end
  if BC == 3; US(1)   = US(2); US(N+2) = US(N+1); end
  

% Time dependent wavefunction UT
  % propagation constant  [1/m]
    k = sqrt(eignV(m,m)) .* (N+1)/L;
  % angular frequency  [rad/s]
    w = v*k;
  % period  [s]
    T = 2*pi/w;
  % frequency  [Hz]
    f = 1/T;
  % time  [s]
    t = linspace(0,2*T,nT);
  % wavelength [m]
    lambda = 2*pi/k;
 % time dependent wavefunction   
    UT = cos(w.*t);   
  
  
 % GRAPHICS  ==========================================================
 figure(1)
  FS = 12;
  pos = [0.05 0.05 0.30 0.20];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  ylim([-1.1 1.1])
    
  xP = x;
  for cc = 1:nT
    yP = US.*UT(cc);
    plot(xP,yP,'b','linewidth',5)
    grid on
    box on
    xlabel('position  x  [ m ]')
    ylabel('wavefunction  U_S  [ a.u.]')
    
    tm1 = 'Mode #  m =  ';
    tm2 = num2str(m,'%3.0f  \n');    
    tm3 = '       \lambda  =  ';
    tm4 = num2str(lambda,'%3.2f m \n');
    tm5 = '       f  =  ';
    tm6 = num2str(f,'%3.1f Hz \n');
    tm = [tm1 tm2 tm3 tm4 tm5 tm6];
    title(tm,'fontsize',FS)
    
    ylim([-1.1 1.1])
    set(gca,'fontsize',FS)
    
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
    
    pause(0.01)
  end

figure(2)
  FS = 12;
  pos = [0.45 0.05 0.30 0.20];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  hold on
  area(x,US,'facecolor',[0.8 0.8 1])
  area(x,-US,'facecolor',[0.8 0.8 1])
  
  xP = x; yP = US;
  plot(xP,yP,'b','linewidth',2)
 
  yP = -US;
  plot(xP,yP,'b','linewidth',2)  
  
  grid on
  set(gca,'fontsize',FS)
  box on
  
  xlabel('position  x  [ m ]')
  ylabel('wavefunction  U_S  [ a.u.]')
  
  tm1 = 'Mode #  m =  ';
  tm2 = num2str(m,'%3.0f  \n');    
  tm3 = '       \lambda  =  ';
  tm4 = num2str(lambda,'%3.2f m \n');
  tm5 = '       f  =  ';
  tm6 = num2str(f,'%3.1f Hz \n');
  tm = [tm1 tm2 tm3 tm4 tm5 tm6];
  title(tm,'fontsize',FS)
 
  
toc
