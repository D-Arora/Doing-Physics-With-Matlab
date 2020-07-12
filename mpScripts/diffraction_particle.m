% diffraction_2slit.m

%171222

% Fraunhofer Diffraction Simulation
% Particle Double Slit Diffraction

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

clear
close all
clc
tic

% INPUTS --------------------------------------------------------

num = 501;                    % No. observation points
num1 = 2 000;
b = 2.0e-4;                   % Slit width
a = 3*b;                    % Slit separation
wL = 400e-9;                  % Wavelength
D = 1.0;                      % Distance to image screen
N = 3;                        % factor for image space range


% SETUP --------------------------------------------------------------
ss = rng;
%rand('state',sum(100*clock))
x1 = -0.004;               % Min position for image space range
x2 = +0.004;               % Max position for image space range
x = linspace(x1,x2,num);   % Data points for image space range

%IRR1 = zeros(num,1);           % Irradiance (intenity)- single slit
%IRR1 = zeros(num,1);           % Irradiance (intenity) - total

k = 2*pi/wL;                  % Wave number
beta = (k*b/(2*D)) .* x;      % Beta parameter
alpha = (k*a/(2*D)) .* x;     % Alpha parameter

% CALCULATIONS --------------------------------------------------

IRR1 = (sin(beta+eps)./(beta+eps)).^2;   % Diffraction - single slit
IRR2 = cos(alpha).^2 .* IRR1;            % Total interference
IRRs = (IRR2.^0.2)' ;     % scaled irradiance for simulated screen pattern
xs = x .* 1e3;               % screen position in mm 
xs1 = x1 * 1e3; xs2 = x2 * 1e3;

xs1 = -4; xs2 = +4;
ys1 = -1; ys2 = +1;

y1 = -1; y2 = +1;

% Graphics ============================================================
figure(1)
set(gcf,'Name','Particle Diffraction');
set(gcf,'NumberTitle','off');
set(gcf,'PaperType','A4');
set(gcf,'Color',[1 1 1]);
fs = 12;

pos = [0.1 0.1 0.33 0.5];
set(gcf,'Units','normalized');
set(gcf,'Position',pos);

% Irradiance - line graph
axes('position',[0.1 0.4 0.8 0.5]);
%plot(xs,IRR1,'k','LineWidth',1);          % Diffraction envelope
hold on
h_axis = area(xs,IRR2);
set(h_axis,'FaceColor',[0 0 0]);
axis([xs1 xs2 0 1.1]);
xlabel('screen position  x  (a.u.)','fontSize',fs);
ylabel('Prob Density (a.u.)','fontSize',fs);
title('Particle Diffraction','fontSize',fs,'fontweight','normal');
set(gca,'fontsize',14)
box on
axes('position',[0.1 0.1 0.8 0.18]);
axis([xs1 xs2 0 1]);
hold on
set(gca,'fontsize',14)

axis off


% CALCULATIONS --------------------------------------------------

% Generate a random number between -1 and +1 for X screen position xP
% Generate a random number between 0 and +1 for Y screen position yP
% Generate a random number between 0 and +1 for screen plotting SP
for cn = 1 : num1
   xP = x1 + (x2-x1)*rand(1,1);
   yP = y1 + (y2-y1)*rand(1,1);
   SP = rand(1,1);
   betaP = (k*b/(2*D)) .* xP;      % Beta parameter
   alphaP = (k*a/(2*D)) .* xP;     % Alpha parameter

   IRRP = (sin(betaP+eps)./(betaP+eps)).^2;   % Diffraction - single slit
   IRRP = cos(alphaP).^2 .* IRRP;            % Total interference

   if SP < IRRP
    plot(xP.*1e3,yP,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','none');
   end
   pause(eps);
end







