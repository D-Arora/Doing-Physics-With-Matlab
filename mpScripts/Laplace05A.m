% Laplace05A.m

% Ian Cooper
% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
% SCRIPTS
%    https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

% SYSTEM: mass m  / spring k / dashpot-damping b.
% Laplace Transform (LT) used to solve the ODE describing the System
%   to give the displacement x(t) of the mass excited by ab input signal f(t).

%   m*x'' + b*x' + k*x = f

% 1   Define the input signal (driving force) f(t).
% 2   LT of input signal --> F(s).
% 3   Define forces acting on System and initial conditions x(0) and x'(0).
% 4   LT of ODE -->
%       m*(s^2*X - s*x(0) - x'(0)) + b*(s*X - x(0)) + k*X = F
%      (m*s^2 + b*s + k)*X - m*(s*x(0) - x'(0)) - b*x(0) - F = 0
%
%  5  Solve equation for X(s)
%  6  Inverse LT --> x(t)

% The symolic solution for the output x(t) is displayed in the
%   Command Window and a summary of results displayed in Figure(2).

% Numeric variables: time tN / displacement xN / velocity vN / acceleration aN

% The Script will need to be amended for different input functions

% 20210303     Matabl R3030b

clear
close all
clc


% System parameters >>>
  m = 1;
  b = 0;
  k = 4;
  x0 = 0;
  v0 = 0;
  w = sqrt(k/m);

  num = 101;
  tMax = 10;
 
  tN = linspace(0,tMax, num);
  
  
syms s t X

% Input function: displacement y and / or driving force f  >>>
%  y = cos(w*t);
%  yN = cos(w*tN);
  
  f = 2;
  fN = 2*ones(num,1);
  
  txt_Input = 'f = -2';
  
  
  
% LT: input  
  F = laplace(f,t,s);


% LT: ODE 
  Z = (m*s^2 + b*s + k)*X - m*(s*x0 - v0) - b*x0 - F;

% Solve for X
  Sol_x = solve(Z, X);

% Inverse LT
  sol_x = ilaplace(Sol_x,s,t)
  pretty(sol_x)
  
% Velocity calculation
  Sol_v = s*Sol_x
  sol_v = ilaplace(Sol_v)  
  

% Extract numbers from symbolic expressions for output
  
  x = subs(sol_x,{t},{tN});
  xN = double(x);
  
  v = subs(sol_v,{t},{tN});
  vN = double(v);
  
% Acceleration calculation  
  aN = -(b/m).*vN - (k/m).* xN + (fN/m)' ;
  
% CALCULATIONS  ==========================================================
% Natural frequency
  wN = sqrt(k/m);
% Gain 
  g = max(xN)/max(fN);
% Natural period
  TN = 2*pi/wN;
% Driving frequency
  T = 2*pi/w;
  
% Period of response - from peaks
  [pks, locs] = findpeaks(xN);
  if length(locs) > 2
     P = ( tN(locs(end)) - tN(locs(2)) )/(length(locs)-2);
  else P = -1;
  end
  

  
% GRAPHICS  ==============================================================
figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.05 0.5 0.30 0.30]);
  set(gcf,'color','w');
  FS = 14; LW = 2;
  box on
  
  plot(tN,xN,'r','linewidth',LW)
  hold on
  
  plot(tN,fN,'b','linewidth',1.2)
  
  grid on
  xlabel('t'); ylabel('x  &  y')
  hL = legend('OUT  x','IN  y');
  set(hL,'orientation','horizontal','location','northoutside','box','off')
  set(gca,'fontsize',FS)

% Results Table: The screen output in Figure(2) can be changed
%                by amending the code as necessary for different inputs.  
figure(2)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.38 0.5 0.25 0.25]);
  set(gcf,'color','w');
  FS = 14; LW = 2;
  box on
  xlim([0 110])
  ylim([0 100])
  
% System 
  hy = 95;
  text(0,hy,'SYSTEM','fontsize',FS);
 
  txt = sprintf('m = %2.2f', m);
  text(30,hy,txt,'fontsize',FS)
  
  txt = sprintf('b = %2.2f', b);
  text(60,hy,txt,'fontsize',FS)
  
  txt = sprintf('k = %2.2f', k);
  text(90,hy,txt,'fontsize',FS);
  
  hy = 80;
  txt = sprintf('\\omega_n = %2.2f', wN);
  text(0,hy-2,txt,'fontsize',FS); 
  
  hy = 80;
  txt = sprintf('T_n = %2.2f', TN);
  text(30,hy-2,txt,'fontsize',FS); 
  
  
  txt = sprintf('x(0) = %2.2f', x0);
  text(60,hy,txt,'fontsize',FS);  
  
  txt = sprintf('v(0) = %2.2f', v0);
  text(90,hy,txt,'fontsize',FS);  

% Input  >>>
  hy = 60;
  text(0,hy,'INPUT','fontsize',FS);
  
  %text(30,hy,'y = cos(\omega t)   f = cos(\omega t))','fontsize',FS);
   text(30,hy,txt_Input,'fontsize',FS);
  
  hy = 42;
  txt = sprintf('\\omega = %2.2f   T = %2.2f', w,T);
 %   text(30,hy,txt,'fontsize',FS);
    

% Response
  hy = 25;
   text(0,hy,'RESPONSE','fontsize',FS);
  
  hy = 10;
  txt = sprintf('Gain  g = %2.2f   P = %2.2f', g,P);
    text(10,hy,txt,'fontsize',FS);  
  
    axis off
    box on
    
