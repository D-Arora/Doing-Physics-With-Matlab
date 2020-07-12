% diffraction_2slit.m

%171222

% Fraunhofer Diffraction: Double Slit

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney

% DOING PHYSICS WITH MATLAB 
%    ../mphome.htm

clear
close all
clc
tic


% INPUTS --------------------------------------------------------
num = 501;                    % No. observation points [501]
b = 1.0e-4;                   % Slit width [2.0e-4]
a = 3*b;                      % Slit separation [3]
wL = 600e-9;                  % Wavelength  [700e-9]
D = 1.0;                      % Distance to image screen [1]
N = 3;                        % factor for image position [3]

% SETUP ---------------------------------------------------------
%xmin = (D/b) * wL;            % Position of 1st min
%x1 = -N * xmin;               % Min position for image space range
%x2 = +N * xmin;               % Max position for image space range
x1 = -10e-3; x2 = 10e-3;
x = linspace(x1,x2,num);      % Data points for image space range

%IRR1 = zeros(num,1);          % Irradiance (intenity)- single slit
%IRR2 = zeros(num,1);          % Irradiance (intenity) - total

k = 2*pi/wL;                  % Wave number
beta = (k*b/(2*D)) .* x;      % Beta parameter
alpha = (k*a/(2*D)) .* x;     % Alpha parameter

% CALCULATIONS --------------------------------------------------
IRR1 = (sin(beta+eps)./(beta+eps)).^2;   % Diffraction - single slit
IRR2 = cos(alpha).^2 .* IRR1;            % Total interference
IRRs = (IRR2.^0.2)' ;     % scaled irradiance for simulated screen pattern
xs = x .* 1e3;                           % screen position in mm 
xs1 = x1 * 1e3; xs2 = x2 * 1e3;

% GRAPHICS ------------------------------------------------------
figure(1)
set(gcf,'Name','Double Slit Diffraction');
set(gcf,'NumberTitle','off');
set(gcf,'PaperType','A4');
graphColor = ColorCode(wL);
set(gcf,'Color',[1 1 1]);
pos = [0.1 0.1 0.38 0.45];
set(gcf,'Units','normalized');
set(gcf,'Position',pos);


% Irradiance - line graph
axes('position',[0.1 0.35 0.78 0.5]);
hplot = plot(xs,IRR1);          % Diffraction envelope
set(hplot,'Color',graphColor,'LineWidth',1);
fs = 12;

hold on
harea = area(xs,IRR2);
set(harea,'EdgeColor',graphColor,'FaceColor',graphColor);
axis([xs1 xs2 0 1.1]);
xlabel('screen position  (mm)','fontSize',fs);
ylabel('Irradiance  ( W.m^{-2} )','fontSize',fs);

mainTitle  = sprintf('Double Slit');
pT1 = '\\lambda = ';
pT2 = num2str(wL,'%3.1e');
pT3 = ' m    b = ';
pT4 = num2str(b,'%3.1e');
pT5 = ' m    a = ';
pT6 = num2str(a,'%3.1e');
pT7 = ' m    a/b = ';  ab = a/b;
pT8 = num2str(ab,'%3.1f');
pT9 = '    D = ';
pT10 = num2str(D,'%3.1f');
pT11 = ' m';

pT = [pT1 pT2 pT3 pT4 pT5 pT6 pT7  pT8 pT9 pT10 pT11];
pT = sprintf(pT);
mainTitle = sprintf('%s\n %s\n ',mainTitle,pT);
h_Title = title(mainTitle);
set(h_Title,'fontweight','normal');
set(gca,'fontsize',14)

% Irradiance -  simulated screen pattern
axes('position',[0.1 0.05 0.8 0.18]);
[xSc, ySc] = meshgrid(xs,[-1,+1]);
zSc = [IRRs,IRRs];
pcolor(xSc', ySc', zSc)
shadingMap = gray(64);
colormap(shadingMap(:,1)*graphColor);
shading interp
axis off








