close all
clear

nfile = 'filament.m';
% calls: simpson1d.m


% *******************************************************************
% DEFINE THE KNOW VARIABLES (SI units, unless stated otherwise)
% *******************************************************************

% Inputs
T = 2400;                      % Temperature of filament
P_total = 55;                  % Total power radiated by filament
num = 201;                     % Number of data points fro wavelength range

% constants
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
wL_peak = b1/T;                % wavelength at peak
wL1 = wL_peak/10;              % min for wavelength range  lambda1
wL2 = 10*wL_peak;              % max for wavelength range  lambda2
wL = linspace(wL1,wL2,num);    % wavelengths

% spectral intensity
P_wL = zeros(num,1);            % array for thermal power per unit wavelength interval

K = (h*c)/(kB*T);              % constant to simply calcuation
P_wL = 1 ./ (wL.^5 .* (exp(K ./ wL)-1));

% Power radiated: area under curve
P = simpson1d(P_wL,wL1,wL2);

% Normalize the power and spectral intensity
N = P/P_total;

P_wL = P_wL/N;

P_check = simpson1d(P_wL,wL1,wL2);  % Check the normalization of the power

%Find visible range
for cn = 1 : num
if wL(cn) > 400e-9
num1 = cn;
break
end; end

for cn = 1 : num
if wL(cn) > 800e-9;
num2 = cn;
break
end; end

% power radiated in visible part of spectrum
P_visible = simpson1d(P_wL,wL(num1),wL(num2));

% efficiency of filament
eff = 100*P_visible/P_total;


% *******************************************************************
% GRAPHICS
% *******************************************************************

figure(1)
Ps = 1e-6;
h_area1 = area(mu*wL,Ps.*P_wL);
hold on
h_area2 = area(mu*wL(num1:num2),Ps.*P_wL(num1:num2));
set(h_area1,'FaceColor',[0 0 0]);
set(h_area2,'FaceColor',[1 1 0]);

tm1 = 'Tungsten filament:  \eta =  ';
tm2 = num2str(eff,3);
tm3 = '  %';
tm = [tm1 tm2 tm3];
hTitle = title(tm);

xlabel('wavelength   \lambda   (\mum)','fontsize',14);
ylabel('Power / d\lambda    (MW.m^{-1})','fontsize',14);
set(gca,'xLim',[0 10]);
set(gca,'fontsize',14);

% *******************************************************************
% SCREEN OUPUT OF RESULTS
% *******************************************************************

clc
disp('    ');
disp(nfile);
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



