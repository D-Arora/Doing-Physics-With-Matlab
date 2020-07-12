function [ ] = tp_fn_Newton(R,N,tMax,T0,Tenv, flagC)

  % Ian Cooper
  % School of Physics, University of Sydney, Sydney, NSW, Australia
  % Email: ian.cooper@sydney.edu.au
  
  % Web:   https://d-arora.github.io/Doing-Physics-With-Matlab/
  % R     Cooling constant [1/min]
  % N     Number of time steps 
  % tMax  Time interval for simulation  [minutes]
  % T0    Initial temperature of system [degC]
  % Tenv  Temperature of surrounding environment  [degc]
  
  % flagC == 1     % Plot numerical solution only for T vs t
  % flagC == 2     % Plot of numerical and analytical solutions for T vs t
  % flagC == 3     % Plot of numerical solution fitted to data for cooling coffee 
  
  close all
  t = linspace(0,tMax,N);      % time
  dt = t(2) - t(1);            % time increment
  T = zeros(1,N);              % ststem temperature as a function of time
  K = dt * R;                  % Constant
  T(1) = T0;                   % Initial temperature of system

% CALCULATIONS ===========================================================

% Numerical computation
   for n = 1:N-1
      T(n+1) = T(n) - K * (T(n) - Tenv);
   end
   Tend = T(end);
   
 % Analytical computation
   TA = Tenv + (T0 - Tenv) .* exp(-R.*t);
 
 % Command Window Output  ===============================================
 disp('   ');
 disp('   ');
 fprintf('Cooling constant               R  = %2.3e   [1/min]  \n',R);
 disp('   ');
 fprintf('Number of time steps           N  = %4.0f  \n',N);
 disp('   ');
 fprintf('Time interval for simulation   tMax  = %4.0f   [min]  \n',tMax);
 disp('   ');
 fprintf('Environmental temperature      Tenv  = %4.2f   [degC]  \n',Tenv);
 disp('   ');
 fprintf('Initial temperature of system  T0  = %4.2f   [degC]  \n',T0);
 disp('   ');
 fprintf('Final temperature of system    Tend  = %4.2f   [degC]  \n',Tend);
 
  
% GRAPHICS ==============================================================

figure(1)
   set(gcf,'units','normalized','position',[0.01 0.2 0.22 0.25]);
   hold on
   xP = t; yP = T;
   plot(xP,yP,'b','lineWidth',2);
   xP = t; yP = Tenv .* ones(1,N);
   plot(xP,yP,'g','lineWidth',1);
   
   xlabel('time  t  [s]')
   ylabel('system   T   [degC]');
   box on
   grid on
   set(gca,'fontsize',14);
   
  
if flagC == 2 
   figure(2)
   set(gcf,'units','normalized','position',[0.25 0.2 0.22 0.25]);
   hold on;
   xP = t; yP = T;
     plot(xP,yP,'b','lineWidth',2);
   xP = t(1:N/10:end); yP = TA(1:N/10:end);
     plot(xP,yP,'ro','lineWidth',0.5);
   xlabel('time  t  [s]')
   ylabel('system   T   [degC]');
   legend('N','A');
   box on
   grid on
   set(gca,'fontsize',14);  
end  % if flagC ==2

if flagC == 3
%    DATA FOR COOLING COFFEE PROBLEM  
     tC = 0:20;    % time recordings
     TC = [90.0 87.0 84.5 82.0 80.0 77.5 75.0 73.0 70.5 68.5 ...
         66.5 65.0 63.5 61.5 59.0 58.5 56.5 55.5 54.0 52.5 51.5];
   
figure(3)
   set(gcf,'units','normalized','position',[0.5 0.2 0.22 0.25]);
   hold on;
   xP = t; yP = T;
     plot(xP,yP,'b','lineWidth',2);
   xP = tC; yP = TC;
     plot(xP,yP,'ro','lineWidth',0.5);
   xlabel('time  t  [s]')
   ylabel('system   T   [degC]');
   %legend('N','A');
   box on
   grid on
   set(gca,'fontsize',14);     
end    % flagC == 3


end

