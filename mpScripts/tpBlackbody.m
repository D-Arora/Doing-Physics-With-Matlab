% tpBlackbody.m

% Blackbody curves for different temperatures

% DOING PHYSICS WITH MATLAB: 
%   ../mphome.htm
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
% Temperatures [K]
  T = [1000 1500 2000 2500]' ;
% Number of data points (must be odd number)
  num = 1501;                     

% Constants
  c = 2.99792458e8;              % speed of light
  h = 6.62608e-34;               % Planck constant
  kB = 1.38066e-23;              % Boltzmann constant
  sigma = 5.6696e-8;             % Stefan constant  
  b1 = 2.898e-3;                 % Wien constant - wavelength
  
  w1 = 2*pi*h*c^2;
  w2 = h*c ./(kB.*T);
  f1 = 2*pi*h/c^2;
  f2 = kB.*T;
  
% *******************************************************************
% CALCULATE THE UNKNOWN VARIABLES (SI uniuts, unless stated otherwise)
% *******************************************************************

% Define the range of values for the wavelength  [m]
  wL_peak = b1./T;                              % wavelength at peak
  wLmin = wL_peak./20;                          % min 
  wLmax = 10.*wL_peak;                          % max 
  wL = linspace(min(wLmin),max(wLmax),num);     % wavelengths
  wL_peak_graph = zeros(length(T),1);
  s = 1e9;                                      % m to nm
  
% Spectral Exitance for wavelength [W.m-2.m-1]   
  R_wL = zeros(length(T),num);  

  for cc = 1 : length(T)
    R_wL(cc,:) = w1 ./ (wL.^5 .* (exp(w2(cc) ./ wL)-1));
  end

 % Power emitted from an area of 1 square meter  A = 1 m^2
   A = 1;
   P = zeros(length(T),1);
   for cc = 1 : length(T)
    P(cc) = A .* simpson1d(R_wL(cc,:),min(wLmin),max(wLmax));
   end
     P = P./min(P);
     
% Wien's Displacement Law   peak wavelength [nm]
 wL_Peak = 1e9.*b1./T;
     
% *******************************************************************
% GRAPHICS
% *******************************************************************

figure(1)

  set(gcf,'units','normalized');
  set(gcf,'position',[0.02 0.45 0.35 0.40]);
  set(gcf,'color','w')
  
  plot(wL.*s,R_wL(1,:),'lineWidth',2,'Color',[1 0 0]);
  hold on
  plot(wL.*s,R_wL(2,:),'lineWidth',2,'Color',[0 0 1]);
  plot(wL.*s,R_wL(3,:),'lineWidth',2,'Color',[1 0 1]);
  plot(wL.*s,R_wL(4,:),'lineWidth',2,'Color',[0 0 0]);
  plot([750 750],[0,13.5e11],'Color',[0.5 0 0]);
  plot([390 390],[0,13.5e11],'Color',[0.3 0 0.3])
 % axis([0 10 0 36])
   xlim([0 5000])
  legend(num2str(T(1)), num2str(T(2)), num2str(T(3)), num2str(T(4)))
  
  grid on
  box on
  title('Thermal Radiation from hot objects');
  xlabel('wavelength   \lambda   [ nm ]');
  ylabel('R_{\lambda}    [ W . m^{-2} . m^{-1} ]');
  text(450,12.5e11,'VIS')
  text(1500,12.5e11,'IR')
  text(100,12.5e11,'UV')
  set(gca,'fontsize',12)

% *******************************************************************
% SCREEN OUPUT OF RESULTS
% *******************************************************************

  disp('    ');
  CWT = table(T,round(P,2),round(wL_Peak,0), 'VariableNames', {'T [K]', 'Prel', 'wL_Peak [nm]'});
  disp(CWT)
  disp('  ')
 

 



