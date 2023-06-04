% mec_osc1.m


close all
clc
clear


% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Forced oscillations
  T0 = 1;
  f0 = 1/T0; w0 = 2*pi*f0;
% max simulation time  [s]
  tMax = 5*T0;
% Initial displacements  [m]
  x0 = 0;
% Initial velocity  v0  [m/s] 
  v0 = 0; 
 

% SETUP ======================================================
  K = [w0];
% tSpan = [0 tMax];
  tSpan = linspace(0,tMax,1201);
% Initial conditions
  u0 = [x0 v0];

% Solve ODE
[t, u] = ode45(@(t,u) EqM(t,u,K), tSpan, u0);

% Extract displacements, velocities
  x = u(:,1);
  v = u(:,2);

% Fourier
  [H, f] = fourier(x',t');
   Hw = H.*hamming(length(H))';

%  RESULTS  ==================================================
%    tF = t(end);       % simulation time
%    vF = sqrt(vx(end)^2 + vy(end)^2);   % final velocity
%    xF = x(end);       % final X displacement
%    yF = y(end);       % final y displacement

% GRAPHICS  ==================================================
figure(1)
   pos = [0.05 0.05 0.27 0.27];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   plot(t,x,'b','linewidth',2)
   grid on
%   hold on
%   txtX = 'x  [ m ]';
%   txtY = 'y  [ m ]'; 
%   txtT = sprintf('t_F = %2.2f   x_F = %2.2f   y_F = %2.2f   v_F = %2.2f \n',tF, xF, yF, vF);
%   plotSetup(txtX,txtY,txtT)

figure(2)
   pos = [0.35 0.05 0.27 0.27];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   plot(f,abs(H),'b','linewidth',2)
   xlim([0 max(f)])
   grid on
% FUNCTIONS  ==================================================

function uDot = EqM(t,u,K)
  w0 = K(1);
  S = 0;
  if t > 1; S = 10; end
  if t > 4; S = 0; end
  uDot = zeros(2,1);
 
  uDot(1) = u(2);
  uDot(2) = -w0^2*u(1) + S ; % + 10*sin(1.5*w0*t) + 10*sin(3*w0*t);  
end

function [H, f] = fourier(h,t)  % Fourier transform H(f)
   nT = length(t); tMin = min(t); tMax = max(t);
   nF = 2999; fMax = 5; fMin = -fMax;
   f = linspace(fMin,fMax,nF);

   H = zeros(1,nF);
   for c = 1:nF
       g = h.* exp(1i*2*pi*f(c)*t);
       H(c) = simpson1d(g,tMin,tMax);
    end
end


