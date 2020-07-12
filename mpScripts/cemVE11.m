% cemVE11.m
% 23 april 2016
% Ian Cooper
% School of Physics, University of Sydney
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% Computation of the electrical potential from the line 
%   integral of the electric field




clear all
close all
clc
tic

% ========================================================================
% INPUTS  
% ========================================================================

% Grid points:  must be an ODD number   
   N = 101;
   
% Number of integration paths
   Np = 3;
% XY path coordinates
   x = [-1.5 -1.5 1 1];
   y = [0 -1 -1 0];
   
% flag = 1: X integration / x inc
% flag = 2: X integration / x dec
% flag = 3: Y integration / y inc
% flag = 4: Y integration / y dec
  flag(1) = 3;
  flag(2) = 1;
  flag(3) = 4;
   
% Dimensions
   minX = -2;  
   maxX =  2;
   minY = -2;
   maxY =  2;

% charge
   Q = 10e-6;
   xC = 0; yC = 0;
% constants
   eps0 = 8.854e-12;
   kC = 1/(4*pi*eps0);
   
% =======================================================================
% CALCULATIONS 
% =======================================================================

V = zeros(Np,1);
% variable pm  plus / minus:  - if x or y increasing / + if x or y
% decreasing

for n = 1 : Np

    if flag(n) == 1
        u = y(n);
        a = x(n); b = x(n+1); pm = -1;
    end
    if flag(n) == 2
        u = y(n);
        a = x(n+1); b = x(n); pm = 1;
    end
   if flag(n) == 3
        u = x(n);
        a = y(n); b = y(n+1); pm = -1;
    end
    if flag(n) == 4
        u = x(n);
        a = y(n+1); b = y(n); pm = 1;
    end 
    
   v = linspace(a,b,N);
   R = sqrt(u.^2 + v.^2);
   R3 = R.^3;
   E = kC .* v .* Q ./R3;
   fn = pm .* E;
   V(n) = simpson1d(fn,a,b); 
end

% Potentials
   V21 = sum(V);        % value of the line integral of E

   R1 = sqrt(x(1)^2 + y(1)^2);
   R2 = sqrt(x(end)^2 + y(end)^2);
   V1 = kC * Q / R1;    % potential at point P1
   V2 = kC * Q / R2;    % potential at point P2
   dV = V2 - V1;        % potential at point P2 w.r.t. point P1

   % display results in Command Window
   disp('  ');
   fprintf('Line integral of E:  V21 = %2.4e  V\n ',V21);
   disp('  ');
   fprintf('Computed potential at point P2 w.r.t point P1  dV = %2.4e  V\n ',dV);
   disp('  ');


% ======================================================================= 
% GRAPHICS 
% =======================================================================

figure(1)   % 11111111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized','position',[0.01 0.52 0.23 0.32]);
   
   xP = xC; yP = yC;
   h = plot(xP,yP,'ro');
   set(h,'markerFacecolor',[1 0 0]);
   hold on
   
   xP = [x(1) x(end)]; yP = [y(1) y(end)];
   h = plot(xP,yP,'bo');
   set(h,'markerFacecolor',[0 0 1]);
   
   xP = x;
   yP = y;
   h = plot(xP,yP,'b');
   set(h,'lineWidth',2);
   
   xlabel('x  [m]'); ylabel('y  [m]');
   
   text(-1.3, -0.3, 'P1','fontsize',12);
   text(1.2 , -0.3, 'P2','fontsize',12);
   text(-0.3, -0.3, 'charge Q','fontsize',12);
   text(-0.2,1.3, '\rightarrow','fontsize',20);
   
   set(gca,'xLim',[minX maxX]);
   set(gca,'yLim',[minY maxY]);
   set(gca,'fontsize',14)

toc