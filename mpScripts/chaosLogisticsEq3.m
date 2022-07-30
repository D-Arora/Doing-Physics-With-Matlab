%  chaosLogisticsEq3.m

% Logistic Difference Equation
% 

% DOING PHYSICS WITH MATLAB
%   https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation
%   https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/???
% Scripts
%   https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%   https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

% IAN COOPER
% email    matlabvisualphysics@gmail.com

% Matlab 2021b     220714

clear
close all
clc

%% CELL 1: PERIOD 1 DYNAMICS
% Control parameter 0 <= r <= 1
  r = 0.3;

% Initial condition   0 <= x0 <= 1
  x0 = 0.95;

% Number of iterations
  N = 5001;

% Logistic Equation
  x = zeros(N,1);
  x(1) = x0;
  
 for c = 1:N-1
     x(c+1) = 4*r*x(c)*(1-x(c));
 end
% Last iterated population
  xEND = x(c+1);
% Fixed (equilibrium) point
  xe = 1-1/(4*r);
  if xe < 0; xe = 0; end

% [1D] iteration of the map f(x) where x(n+1) = 4*pi*r*(1-x(n))
  Nx = 501;
  xn = linspace(0,1,Nx);
  fn = 4.*r.*xn.*(1-xn);

  dx = xn(2)-xn(1);
  gradF = gradient(fn,dx);
   

figure(1)
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.05 0.05 0.25 0.30])
   set(gcf,'color','w');
   FS = 14;
 
% [1D] map  f(x) vs x
  xP = xn;
  yP = fn;
  plot(xP,yP,'b-','linewidth',2)

  hold on
%  y = x
   yP = xP;
   plot(xP,yP,'r-','linewidth',2)

   yP = abs(gradF);
   plot(xP,yP,'m-','linewidth',2)
   xlabel('x')

   xP = [xe xe]; yP = [0 abs(4*r - 8*r*xe)];
   plot(xP,yP,'m','linewidth',0.5)
   
   xP = [0 xe]; yP = [abs(4*r - 8*r*xe) abs(4*r - 8*r*xe)];
   plot(xP,yP,'m','linewidth',0.5)

   ylabel('f(x)')
   txt = sprintf('r = %2.3f  x_e = %2.4f   x_{END} = %2.4f   \n',r, xe, xEND);
   title(txt)
   xlim([0 1]);
   ylim([0 1.2]); 
   set(gca,'fontsize',FS)
   grid on
   box on

%% CELL 2: PERIOD 1 DYNAMICS
%  FIXED POINT xe AS A FUNCTION OF GROWTH PARAMETER r  ===============

  rMin = 0.25; rMax = 0.75;
  r = linspace(rMin,rMax, 201);
  xe = 1 - 1./(4.*r);

figure(2)
  set(gcf,'units','normalized');
  set(gcf,'Position', [0.31 0.05 0.25 0.30])
  set(gcf,'color','w');
  FS = 14;
  
  xP = r; yP = xe;
   plot(xP,yP,'b-','linewidth',2)
   hold on
   
   xP = [0.25 0.25]; yP = [0 1];
    plot(xP,yP,'k-','linewidth',1)
   xP = [0.75 0.75]; yP = [0 1];
    plot(xP,yP,'k-','linewidth',1)

  xP = [0 0.75]; yP = [xe(end) xe(end)];
    plot(xP,yP,'r-','linewidth',1)

   xlabel('r')
   ylabel('x_e')
   txt = sprintf('max(x_e) = %2.4f   \n',xe(end));
   text(0.4, 0.7,txt,'Color','r','FontSize',FS)

   txt = 'non-trival fixed point x_e';
   title(txt,'FontWeight','normal')
   xticks(0.2:0.1:0.8)
   xlim([0.2 0.8]); ylim([0 1]);
   set(gca,'fontsize',FS)
   grid on
   box on



%%  CELL 3: PERIOD 2 DYNAMICS 

% Growth parameter r  0 < r < 1
  r = 0.80;  

% Period
  period = 2;
 
% Initial population  
  x0 = 0.7;

 if period == 2
    if r > 0.855; r = 0.85; end
 end

x = linspace(0.01,0.99,1201);
dx = x(2) - x(1);
%x1 = 4*r.*x.*(1-x);
x1 = iterate(x,r);
x2 = iterate(x1,r);
x3 = iterate(x2,r);
x4 = iterate(x3,r);

gradF = gradient(x2,dx);
% for c = 2:N
%     f(c) = iterate(f(c-1),r);
% end


% Intersection points x and f(x)   0 < r < 0.855
if period == 2
  interSections = find(abs(x2-x) <= 0.001);
  xIS = x(interSections);
  disp('intersection points')
  fprintf('xIS = %2.3f   ',xIS)
  disp('  ')
  
end  

figure(3)
 set(gcf,'units','normalized');
  set(gcf,'Position', [0.57 0.05 0.25 0.30])
  set(gcf,'color','w');
  FS = 14;

  xP = x;
  if period == 1; yP = x1; end
  if period == 2; yP = x2; end
  if period == 3; yP = x3; end
  if period == 4; yP = x4; end

  hold on

  plot(xP,yP,'b','linewidth',2)

  plot(x,x,'r-','linewidth',2)

  plot(x,abs(gradF),'k')
  
  ylim([0 1.2])
  grid on
  box on
  xlabel('x')
  ylabel('y')
  txt = sprintf('r = %2.3f  \n',r);
  title(txt)
  set(gca,'fontsize',FS)

% ====================================================================
function f = iterate(x,r)

f = (4*r).*x.*(1-x);
 
end


