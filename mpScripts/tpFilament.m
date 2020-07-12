% tpFilamant.m

% Blackbody: Hot tungsten filament 

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/tpBlackbody.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper
% email: matlabvisualphysics@gmail.com
% 191206

% Calls: simpson1d.m  ColorCode.m


close all
clear
clc


% *******************************************************************
% DEFINE THE KNOW VARIABLES (SI units, unless stated otherwise)
% *******************************************************************

% Inputs
  T = 2400;                      % Temperature of filament  [K]
  P_total = 55;                  % Total power radiated by filament [W]
  num = 801;                     % Number of data points fro wavelength range

% constants
  c = 2.99792458e8;              % speed of light [m/s]
  h = 6.62608e-34;               % Planck constant [J.s]
  kB = 1.38066e-23;              % Boltzmann constant 
  sigma = 5.6696e-8;             % Stefan constant  
  b1 = 2.898e-3;                 % Wien constant - wavelength
  mu = 1e6;                      % convert m into micrometers


% *******************************************************************
% CALCULATE THE UNKNOWN VARIABLES (SI uniuts, unless stated otherwise)
% *******************************************************************

% Define the range of values for the wavelength 
  wL_peak = b1/T;                % wavelength at peak  
  wL1 = wL_peak/10;              % min for wavelength range  lambda1
  wL2 = 10*wL_peak;              % max for wavelength range  lambda2
  wL = linspace(wL1,wL2,num);    % wavelengths [m]

% spectral exitance  * area   [W.m-1]
  R_wL = zeros(num,1);            % array for thermal power per unit wavelength interval

  K = (h*c)/(kB*T);              % constant to simply calcuation
  R_wL = 1 ./ (wL.^5 .* (exp(K ./ wL)-1));

% Power radiated: area under curve  [W]
  P = simpson1d(R_wL,wL1,wL2);

% Normalize the power 
  N = P/P_total;

  R_wL = R_wL/N;   % [W/m]

  P_check = simpson1d(R_wL,wL1,wL2);  % Check the normalization of the power

%Find visible range
for cn = 1 : num
  if wL(cn) > 400e-9
    num1 = cn;
    break
  end
end

for cn = 1 : num
  if wL(cn) > 800e-9
     num2 = cn;
     break
  end
end

% Power radiated in visible part of spectrum
  P_visible = simpson1d(R_wL,wL(num1),wL(num2));

% Efficiency of filament
  eff = 100*P_visible/P_total;


% *******************************************************************
% GRAPHICS
% *******************************************************************

figure(1)
  set(gcf,'units','normalized');
 set(gcf,'position',[0.02 0.05 0.3 0.35]);
 set(gcf,'color','w');
 
  Ps = 1e-6;
  h_area1 = area(mu*wL,Ps.*R_wL);
  hold on
  h_area2 = area(mu*wL(num1:num2),Ps.*R_wL(num1:num2));
  set(h_area1,'FaceColor',[0 0 0]);
  set(h_area2,'FaceColor',[1 1 0]);

  tm1 = 'Tungsten filament:  \eta =  ';
  tm2 = num2str(eff,2);
  tm3 = '  %';
  tm = [tm1 tm2 tm3];
  hTitle = title(tm);

  xlabel('wavelength   \lambda   [ \mum ]','fontsize',14);
  ylabel('Power / d\lambda    [ MW.m^{-1} ]','fontsize',14);
  set(gca,'xLim',[0 5]);
  set(gca,'fontsize',14);
  
  

% *******************************************************************
% SCREEN OUPUT OF RESULTS
% *******************************************************************

disp('    ');
fprintf('   wavelength at peak  = %0.2e  m  \n',wL_peak); 
fprintf('   wavelength at peak  = %0.2f  um  \n',mu*wL_peak); 
disp('    ');
fprintf('   P_total    = %0.1f  W  \n',P_total);
disp('    ');
fprintf('   P_visible  = %0.1f  W  \n',P_visible); 
disp('    ')
fprintf('   efficiency (percentage) =   %0.1f  \n',eff); 
disp('    ');
disp('Check normalization');
fprintf('   P_check    = %0.1f  W  \n',P_check);



