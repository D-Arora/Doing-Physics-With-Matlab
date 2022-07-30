%  chaosLogisticsEq03.m

% Logistic Difference Equation
% Iterations of map f(x) = x(n+1) = 4*r*x(n)*(1-x(n))
% PERIOD DOUBLING

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


% INPUTS ==============================================================
% Period  p = 1, 2, ..., 8
  p = 6;

% Control parameter 0 <= r <= 1
  r = 0.960;

% Number of iterations
  N = 100001;

% SETUP  =============================================================  
% Population
x = linspace(0,1,N);
  dx = x(2) - x(1);

% Logistic Difference Equations for period doubling
  f1 = iterate(x,r);
  f2 = iterate(f1,r);
  f3 = iterate(f2,r);
  f4 = iterate(f3,r);
  f5 = iterate(f4,r);
  f6 = iterate(f5,r);
  f7 = iterate(f6,r);
  f8 = iterate(f7,r);

  if p == 1; f = f1; end
  if p == 2; f = f2; end
  if p == 3; f = f3; end
  if p == 4; f = f4; end
  if p == 5; f = f5; end
  if p == 6; f = f6; end
  if p == 7; f = f7; end
  if p == 8; f = f8; end

% Gradient df/dx
  gradF = gradient(f,dx);

% FIXED POINTS: Intersection points x and f(x)   
  interSections = find(abs(f-x) <= 0.00001);
%  xIS = x(interSections);
 
  n = 0; xe(1) = 0;
  for c = 1 : length(interSections)
   if abs(gradF(interSections(c))) < 1
     n = n+1;
     xe(n) = x(interSections(c));
   end
  end
 gradF_xe = zeros(length(c),1);
 for c = 1: length(xe)
     z = find(x > xe(c),1);
     gradF_xe(c) = gradF(z);
 end
 
 disp('Intersection points')
 for c = 1: length(xe)
     fprintf('xe = %2.4f   df(xe)/dx = %2.4f   \n',xe(c), gradF_xe(c))
 end
 disp('  ')


 
% GRAPHICS  =======================================================  
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'Position', [0.05 0.55 0.45 0.30])
   set(gcf,'color','w');
   FS = 14;
 
% Plot F(x) vs x  
  xP = x; yP = f;
  plot(xP,yP,'b-','linewidth',2)

  hold on
%  Plot straight line f(x) = x
   yP = xP;
   plot(xP,yP,'r-','linewidth',2)
% Plot gradient F(x): set max to 1.1
   yP = abs(gradF);
   yP(yP>1.1) = 1.001;
   plot(xP,yP,'m-','linewidth',0.25)

   xlabel('x'); ylabel('f(x)')
   txt = sprintf('r = %2.3f\n',r);
   title(txt)
   xlim([0 1]); ylim([0 1.001]); 
   set(gca,'fontsize',FS)
   grid on; box on

% PERIOD 1 DYNAMICS  ================================================   
% FIXED POINT xe AS A FUNCTION OF GROWTH PARAMETER r  
if p == 1
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
   grid on; box on
end


% %%  CELL 3: PERIOD 2 DYNAMICS 
% 
% % Growth parameter r  0 < r < 1
%   r = 0.80;  
% 
% % Period
%   period = 2;
%  
% % Initial population  
%   x0 = 0.7;
% 
%  if period == 2
%     if r > 0.855; r = 0.85; end
%  end
% 
% x = linspace(0.01,0.99,1201);
% dx = x(2) - x(1);
% %x1 = 4*r.*x.*(1-x);
% x1 = iterate(x,r);
% x2 = iterate(x1,r);
% x3 = iterate(x2,r);
% x4 = iterate(x3,r);
% 
% gradF = gradient(x2,dx);
% % for c = 2:N
% %     f(c) = iterate(f(c-1),r);
% % end
% 
% 
% % Intersection points x and f(x)   0 < r < 0.855
% if period == 2
%   interSections = find(abs(x2-x) <= 0.001);
%   xIS = x(interSections);
%   disp('intersection points')
%   fprintf('xIS = %2.3f   ',xIS)
%   disp('  ')
%   
% end  
% 
% figure(3)
%  set(gcf,'units','normalized');
%   set(gcf,'Position', [0.57 0.05 0.25 0.30])
%   set(gcf,'color','w');
%   FS = 14;
% 
%   xP = x;
% %   if period == 1; yP = x1; end
% %   if period == 2; yP = x2; end
% %   if period == 3; yP = x3; end
% %   if period == 4; yP = x4; end
% 
%   hold on
% 
%   plot(xP,yP,'b','linewidth',2)
% 
%   plot(x,x,'r-','linewidth',2)
% 
%   plot(x,abs(gradF),'k')
%   
%   ylim([0 1.2])
%   grid on
%   box on
%   xlabel('x')
%   ylabel('y')
%   txt = sprintf('r = %2.3f  \n',r);
%   title(txt)
%   set(gca,'fontsize',FS)

% ====================================================================
function f = iterate(x,r)

f = (4*r).*x.*(1-x);
 
end


