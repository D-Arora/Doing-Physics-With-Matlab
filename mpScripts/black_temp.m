% black_temp.m

% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/tp_blackbody.pdf
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% Ian Cooper  MatlabVisualPhysics@gmail.com
% 191205

% calls: simpson1d.m
% Plots - spectral intensity curves, 4 temperatures
% Calculations - Peak wavelengths, Total power radiated


close all
clear
clc

% *******************************************************************
% DEFINE THE KNOW VARIABLES (SI units, unless stated otherwise)
% *******************************************************************

% Inputs
  T = [1000 1500 2000 2500]     % Temperatures [K]
  P_total = 100;                  % Total power radiated by filament
  num = 201;                     % Number of data points fro wavelength range

% Constants
  c = 2.99792458e8;              % speed of light
  h = 6.62608e-34;               % Planck constant
  kB = 1.38066e-23;              % Boltzmann constant
  sigma = 5.6696e-8;             % Stefan constant  
  b1 = 2.898e-3;                 % Wien constant - wavelength
  b2 = 2.82*kB/h;                % Wien constant - frequency
  mu = 1e6;                      % convert m into micrometers


% *******************************************************************
% CALCULATE THE UNKNOWN VARIABLES (SI uniuts, unless stated otherwise)
% *******************************************************************

% Define the range of values for the wavelength 
  wL_peak = b1./T;                         % wavelength at peak
  wL1 = wL_peak./20;                       % min for wavelength range  lambda1
  wL2 = 10.*wL_peak;                       % max for wavelength range  lambda2
  wL = linspace(min(wL1),max(wL2),num);    % wavelengths

  P_wL = zeros(length(T),num);  % array for thermal power per unit wavelength interval

  wL_peak_graph = zeros(length(T),1);
  N = zeros(length(T),1);
  P = zeros(length(T),1);

% Spectral Intensity
for cc = 1 : length(T)
 K = (h*c)/(kB*T(cc));                 % constant to simply calcuation
 P_wL(cc,:) = 1 ./ (wL.^5 .* (exp(K ./ wL)-1));
end

for cc = 1 : length(T)
  N(cc) = max(P_wL(cc,:));
  P_wL(cc,:) = P_wL(cc,:)/N(1);
  wL_peak_graph(cc) = wL(P_wL(cc,:) == max(P_wL(cc,:))); 

% Power radiated: area under curve
  funct = P_wL(cc,:);
  P(cc) = simpson1d(funct,min(wL1),max(wL2));
end
  P = P.*P_total./max(P);
  

% *******************************************************************
% GRAPHICS
% *******************************************************************

figure(1)

  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.45 0.35 0.40]);
  set(gcf,'color','w')
  
  plot(wL*1e6,P_wL(1,:),'lineWidth',2,'Color',[1 0 0]);
  hold on
  plot(wL*1e6,P_wL(2,:),'lineWidth',2,'Color',[0 0 1]);
  plot(wL*1e6,P_wL(3,:),'lineWidth',2,'Color',[1 0 1]);
  plot(wL*1e6,P_wL(4,:),'lineWidth',2,'Color',[0 0 0]);
  plot([0.75 0.75],[0,36],'Color',[0.5 0 0]);
  plot([0.39 0.39],[0,36],'Color',[0.3 0 0.3])
 % axis([0 10 0 36])
  xlim([0 10])
  legend(num2str(T(1)), num2str(T(2)), num2str(T(3)), num2str(T(4)))
  
  grid on
  box on
  title('Thermal Radiation from hot objects');
  xlabel('wavelength   \lambda   (\mum)');
  ylabel('Power / d\lambda    (W.m^{-1})');
  text(0.4,32,'VIS')
  text(3,32,'IR')
  set(gca,'fontsize',12)

% *******************************************************************
% SCREEN OUPUT OF RESULTS
% *******************************************************************

  disp('    ');
  disp('Temperatures')
  fprintf('  %0.0f  K  \n',T);
  disp('    ');
  disp('Relative Peak wavelengths ')
  fprintf('  %0.1f   \n',wL_peak_graph/wL_peak_graph(1));

  disp('    ');
  disp('Relative total power radiated (area under curves)')
  fprintf('  %0.1f   \n',P/P(1));




