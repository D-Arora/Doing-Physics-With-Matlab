% chaos10.m

% Dynamics of Linear and Nonlinear Systems
% PHASE PLANE ANALYSIS: LInear Systems

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180813 / Matlab version R2018a

% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Reference page for documentation and notes
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/chaos10.pdf

clear 
close all  % comment/uncomment to hold on plot for extra initial conditions
clc


% INPUTS [default values] ===========================================

% D.E. coefficients  K matrix  [2 1 0; 2 -1 0]
  k11 = -40; k12 = -160; k13 = 400;
  k21 = 12.5; k22 = -12.5; k23 = 0;

% Initial conditions t = 0 : range for x1(1)  [-10:10]
   x1I = 0;
   
% Initial conditions: t = 0  x2(1) values  [9]
   x2I = 0;  

% Time domain: time interval tMax [15] / number of time steps nT [2000]
   tMax = 0.5;
   nT = 2000;

% Phase space setup: dimensions [-10 to + 10] / number of vectors nX*nX [16] 
   L = 10;
   x1Min = -L; 
   x1Max = L;
   x2Min = -L; 
   x2Max = L;
   nX = 16;      

     
% CALCULATIONS =====================================================

% K matrix
   K = [k11 k12;k21 k22];
 
% Equilibrium (critical) point 
   xC = K\[-k13;-k23];
    
% (x,y) coordinates for a set of Initial Conditions
   LenXI = length(x1I);
   x1I(LenXI+1:2*LenXI) = x1I;
   x2I = x2I.*ones(1,LenXI);
   x2I(LenXI+1:2*LenXI) = -x2I;

% Eigenfunctions a and eigenvalues b
   [a, b] = eig(K);
   

% COMMAND WINDOW OUTPUT  ==============================================
  disp('D.E. coefficients k11 k12 k13 / k21 k23 k23')
  fprintf('   %2.2f   %2.2f   %2.2f ',k11,k12,k13);
  disp('  ')
  fprintf('   %2.2f   %2.2f   %2.2f ',k21,k22,k23);
  disp('  ')
  disp('  ')
  disp('Eigenvalues b = ')
  disp(b)
  disp('Eigenfunction a = ')
  disp(a)
  
  
% SOLUTIONS D.E.  =====================================================  
figure(1)   % phase space portrait
   FS = 12;
   pos = [0.02 0.05 0.25 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
for c = 1:2*LenXI 
% Initial conditions    
   xI = [x1I(c);x2I(c)];

% C coefficients
   C = a\(xI-xC);

% c coefficients
   cc = zeros(2,2);
   cc(:,1) = a(:,1)*C(1);
   cc(:,2) = a(:,2)*C(2);

% Time domain 
   t = linspace(0,tMax,nT);
  
% Phase space state variables: Solutions  
   x1 = cc(1,1).*exp(b(1,1)*t) + cc(1,2).*exp(b(2,2)*t);
   x2 = cc(2,1).*exp(b(1,1)*t) + cc(2,2).*exp(b(2,2)*t);

% Translations  
   x1 = x1 + xC(1);
   x2 = x2 + xC(2);
  
% Phase plane grid and unit vectors
   x = linspace(x1Min,x1Max,nX);
   y = linspace(x2Min,x2Max,nX);
   [xx, yy] = meshgrid(x,y);
   f = K(1,1).*xx + K(1,2).*yy + k13;
   g = K(2,1).*xx + K(2,2).*yy + k23;
  
   fs = f./sqrt(f.^2 + g.^2);    % unit vectors
   gs = g./sqrt(f.^2  +g.^2);
  
% x1 and x2 nullclines: slopes m 
   mx = -(K(1,1))/K(1,2);
   my = -(K(2,1))/K(2,2);
   xN = mx .* x - k13/k12;
   yN = my .* x - k23/k22;
   


   
% GRAPHICS ============================================================   

% VECTOR FIELD
   hq = quiver(xx,yy,fs,gs);
   xlim([x1Min x1Max])
   ylim([x2Min x2Max])
   
   set(hq,'color',[0.1 0.1 0.1],'AutoScaleFactor',0.6);
   set(gca,'fontsize',FS)
   
   hold on
   
% NULLCLINES
   xP = x; yP = xN;
   if K(1,2) == 0
     xP = [xC(1) xC(1)];  yP = [x2Min x2Max]; 
   end    
   plot(xP,yP,'r');
   
   xP = x; yP = yN;
   if K(2,2) == 0
     xP = [xC(1) xC(1)];  yP = [x2Min x2Max]; 
   end    
   plot(xP,yP,'m');
   
% TRAJECTORY   
   plot(real(x1),real(x2),'b','linewidth',2);
   
% INITIAL CONDITIONS   
   Hplot = plot(real(x1(1)),real(x2(1)),'o');
   set(Hplot,'markersize',8,'markerFaceColor',[0.5 0.5 1],'markerEdgeColor',[0.5 0.5 1])
   
% Graph setup   
   axis square
   xlabel('x_1')
   ylabel('x_2')
   box on
   tm1 = 'b_{11} = ';
   tm2 = num2str(b(1,1),'%3.2f ');
   tm3 = '   b_{22} = ';
   tm4 = num2str(b(2,2),'%3.2f');
   tm = [tm1 tm2 tm3 tm4];
   ht = title(tm,'fontweight','normal');
   set(ht,'Fontsize',14);
end

figure(2)  % state variables as functions of time
   FS = 12;
   pos = [0.32 0.05 0.25 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   
   xP = t; yP = real(x1);
   plot(xP,yP,'r','linewidth',2);
   hold on
   yP = real(x2);
   plot(xP,yP,'b','linewidth',2)
   grid on
   legend('x_1','x_2');
   xlabel('t'); ylabel('x');
   tm1 = '  t = 0:  x_1(0) = ';
   tm2 = num2str(real(x1(1)),'%3.2f ');
   tm3 = '   x_2(0) = ';
   tm4 = num2str(real(x2(1)),'%3.2f');
   tm = [tm1 tm2 tm3 tm4];
   title(tm,'fontweight','normal')
   set(gca,'fontsize',FS)
   
   
   