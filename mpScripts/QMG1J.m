% QMG1H.m

% DOING PHYSICS WITH MATLAB
% https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
% https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/QMG102.htm
% IAN COOPER
% matlabvisualphysics@gmail.com
% 230415   Matlab R2021b


% GRIFFITH  problem 1.13
clear; close all; clc
% SETUP
  A = 1;                    % amplitude
  k = 2;                    % spring constant
  m = 2;                    % mass
  w = sqrt(k/m);            % angular frequency
  T = 2*pi/w;               % period
% SIMULATION: randon times
  N = 1e5+1;                  % Number of random times
  t = (T/1).*rand(1,N);
  x = A.*cos(w*t);
  v = -w.*A.*sin(w*t);
  p = m.*v;
  xmean = mean(x);
  xstd = std(x);
  pmean = mean(p);
  pstd = std(p);
% PROBABILITY and EXPECTATION VALUES
  xavg = simpson1d(x,0,T)/T;
  x2avg = simpson1d(x.^2,0,T)/T;
  sigmaX = sqrt(x2avg - xavg^2);
  pavg = simpson1d(p,0,T)/T;
  p2avg = simpson1d(p.^2,0,T)/T;
  sigmaP = sqrt(p2avg - pavg^2);
  rho = 1./(pi.*sqrt(A^2 - x.^2));
  syms xs
  checkRHO = int(1/(pi*sqrt(1-xs^2)),[-1 1]);
  
% OUTPUT 
    disp('S simulation results / T probability & expectation results')
    fprintf('S number of random times, N = %2.4e  \n', N)
    fprintf('S average position, <x> = %2.4f  \n', xmean)
    fprintf('T average position, <x> = %2.4f  \n', xavg)
    fprintf('S standard deviation position, sigmaX = %2.4f  \n', xstd)
    fprintf('T standard deviation position, sigmaX = %2.4f  \n', sigmaX)
    fprintf('S average momentum, <p> = %2.4f  \n', pmean)
    fprintf('T average momentum, <p> = %2.4f  \n', pavg)
    fprintf('S standard deviation momentum,  sigmaP = %2.4f  \n', pstd)
    fprintf('T standard deviation momentum,  sigmaP = %2.4f  \n', sigmaP)
    fprintf('T CHECK rho normalization = %2.4f  \n', checkRHO)

% GRAPHICS
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.05 0.25 0.60]);
  set(gcf,'color','w');
  FS = 14;
subplot(4,1,1)
  xP = t; yP = x;
  plot(xP,yP,'b.','LineWidth',2)
  %xticks([L1:1:L2])
  grid on
  xlabel('t')
  ylabel('x')
  title('position','FontWeight','normal')
  set(gca,'FontSize',FS)
subplot(4,1,2)
  xP = t; yP = p;
  plot(xP,yP,'b.','LineWidth',2)
  %xticks([L1:1:L2])
  grid on
  xlabel('t')
  ylabel('p')
  set(gca,'FontSize',FS)
  title('momentum','FontWeight','normal')
subplot(4,1,3)
  R = linspace(-1,1,998);
  dR = R(2)-R(1);
  h = histogram(x,'BinEdges',R);
  c = get(h,'Values');
  pD = c./(N*dR);
  check = simpson1d(pD,-1,1);
  ylim([0 1000])
  grid on
  xlabel('x')
  ylabel('\rho(x)')
  set(gca,'FontSize',FS) 
  title('normalized porbability density','FontWeight','normal')
  fprintf('S Check normalization prob. density = %2.4f  \n', check)
subplot(4,1,4)
  plot(x,rho,'.')
  ylim([0 2])
  grid on
  xlabel('x')
  ylabel('\rho(x)')
  title('\rho(x) = 1/(\pi\surd(1-x^2))','FontWeight','normal') 
  set(gca,'FontSize',FS)  
