% op_bessel1.m
% Bessel function of the first kind J1
% 30 aug 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% Plot of the Bessel function of the first kind
% Plots of the Fraunhoffer diffraction pattern for a uniformly illuminated
%    circular aperture
% Calls turningPoint.m  to find indices for the zeros, max, min of a function


% ************************************************************************
clear all
close all
clc

N = 99999;                    % no. of partitions
vMin = 0; vMax = 50;         % optical radial coordinate
v = linspace(vMin, vMax, N);

% Bessle Function ---------------------------------------------------------
J1 = besselj(1,v);           % Bessel function of the first kind
y = max(J1) .* sin(v);       % Sine function

% Turning Points for Bessel Function   min  max  zeros
    xData = v; yData = J1;  
    [indexMin, indexMax] = turningPoints(xData, yData);
    disp('J1   BESSEL FUNCTION OF THE FIRST KIND  ')
    disp('MINs: Radial coordinate  / J1 value');
    for c = indexMin
       fprintf('     %3.3f   %3.3f  \n' ,v(c), J1(c));
    end
    disp('  ')
    disp('MAXs: Radial coordinate  / J1 value');
    for c = indexMax
       fprintf('     %3.3f   %3.3f  \n' ,v(c), J1(c));
    end
    disp('  ')
    yData = J1.^2;  
    [indexMin, indexMax] = turningPoints(xData, yData);
    disp('ZEROs: Radial coordinate    J1 = 0');
    for c = indexMin
       fprintf('     %3.3f   \n' ,v(c));
    end
    disp('  ')
   
    
% ========================================================================    
% Fraunhoffer diffraction from a circular aperture -----------------------
IRR = (J1 ./ v).^2;          % irradiance      
IRR = IRR ./ max(IRR);
IRRs = 10 .* log10(IRR);

% Turning Points for INTENSITY  max  zeros
    xData = v; yData = IRR; 
    [indexMin, indexMax] = turningPoints(xData, yData);
    disp('IRR Fraunhofer Diffraction   MAX   ZEROS  ')
    disp('MAXs: Radial coordinate  / IRR value');
    for c = indexMax
       fprintf('     %3.3f   %3.5f  \n' ,v(c), IRR(c));
    end
    disp('  ')
    disp('ZEROs: Radial coordinate    IRR = 0');
    for c = indexMin
       fprintf('     %3.3f   \n' ,v(c));
    end
    disp('  ')


% ========================================================================

% Graphics ---------------------------------------------------------------
figure(1)                    % Bessel and sine fucntions
fs = 12;
tx = 'v';
ty = 'J_1';
set(gca,'fontsize',fs);
set(gcf,'Units','normalized');
set(gcf,'Position',[0.2 0.2 0.3 0.3]);
plot(v,J1,'linewidth',2)
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
grid on
hold on
plot(v,y,'m','linewidth',1)
legend('J_1','sine');

figure(2)                    % Irradiance
fs = 12;
tx = 'v';
ty = 'IRR (dB)';
set(gca,'fontsize',fs);
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.2 0.3 0.4]);
subplot(2,1,1);
plot(v,IRR,'linewidth',2)
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
grid on
subplot(2,1,2);
plot(v,IRRs,'linewidth',2);
xlabel(tx,'fontsize',fs);
ylabel(ty,'fontsize',fs);
grid on

