% spCalculations.m

%%  mod52A.htm
% radius
   r = 0.74;
% Force magnitudes
  F = [45 54 75];
% Angles
  theta = [90 50 75];
  
% Torques  
  tau = r .* F .* sind(theta);
  
%% mod53A
mA = 0.5;
mB = 0.01
R = 0.75*2
G = 6.67408e-11

F1 = G*mA*mB/R^2

aA = F1/mA
aB = F1/mB

  
  
  