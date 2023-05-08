% QMG1B.m

% Run CELLS #! and #2 separately

%%   CELL #1
clear;close all;clc

% SETUP
  h = 1;
  g = 9.8;
  nT = 1e6;

 T = sqrt(2/g);      % time of flight
  t = T.*rand(nT,1); % time
 x = 0.5*g.*t.^2;    % height measurement
 xavg = mean(x);     % <x>    expectation value
 xstd = std(x);      % sigma  standard deviation 

% Probability <x> - sigma <= x <= <x> + sigma
  ind1 = find(sort(x) > xavg-xstd,1);
  ind2 = find(sort(x) > xavg+xstd,1);
  prob1 = (ind2-ind1)/nT;
% Probability x < <x> - sigma and x > <x> + sigma
  prob2 = 1 - prob1;

% OUTPUT
  fprintf('time of flight, T = %2.4f  \n', T)
  fprintf('expectation value, <x> = %2.4f  \n', xavg)
  fprintf('standard deviation of x, std = %2.4f  \n', xstd)
  fprintf('Probability <x> - sigma <= x <= <x> + sigma, prob1 = %2.4f  \n', prob1)
   fprintf('Probability x < <x> - sigma and x > <x> + sigma , prob2 = %2.4f  \n', prob2)

 % GRAPHICS 
   figure(1)
   pos = [0.52 0.05 0.22 0.22];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   plot(t,x,'.')
   xlabel('t')
   ylabel('x')
   set(gca,'FontSize',14)
   grid on


%%  CELL #2
clear;close all;clc
syms x h 

% probability density
  probD = 1/(2*sqrt(h*x));
% probabilty symbolic  
  Ps = int(probD,x);
% Probability  
  P = int(probD,x,[0 h]);
  fprintf('Probability 0<= x <= 1), P = %2.4f \n',P)
% Expection value for position <x>
  xProbD = x*probD;
  xavg = int(xProbD,x,[0 h]);
  disp(' ')
  displayFormula(" xAVG = xavg")
% Standard deviation  h = 1
  x2avg =  int(x^2*probD,[0 h]);
  displayFormula("x2AVG = x2avg")
  xavg2 = xavg^2;
  displayFormula("xAVG2 = xavg2")
  sigma2 = x2avg - xavg2;
  displayFormula("variance = sigma2")
  sigma = sqrt(sigma2);
  displayFormula(" STD = sigma")
  fprintf('STD/h = %2.4f \n',subs(sigma,1))
  disp('  ')
% Prob within 1 STD of mean
  %L1 = subs(xavg -xSTD,1); L2 = subs(xavg + xSTD,1);
  L1 = subs(xavg - sigma,1); L2 = subs(xavg + sigma,1);
  Psigma = (int(probD,[L1 L2]));
  disp('Prob <x> - sigma <= x <= <x> + sigma')
  fprintf('  Psigma = %2.4f \n',subs(Psigma,1))
  disp('Prob x < <x> - sigma or x > <x> + sigma')
  fprintf('  1 - Psigma = %2.4f \n',1 - subs(Psigma,1))


  
