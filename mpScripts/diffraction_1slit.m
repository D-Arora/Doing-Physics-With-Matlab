% diffraction_2slit.m

%171222

% Fraunhofer Diffraction: Single Slit

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
num = 501;                    % No. observation points
b = 2.0e-4;                   % Slit width
wL = 700e-9;                  % Wavelength
D = 1.0;                      % Distance to image screen
N = 4;                        % factor for image space range

% SETUP ---------------------------------------------------------
xmin = (D/b) * wL;            % Position of 1st min
x1 = -N * xmin;               % Min position for image space range
x2 = +N * xmin;               % Max position for image space range
x = linspace(x1,x2,num);      % Data points for image space range

IRR = zeros(num,1);           % Irradiance (intenity)

k = 2*pi/wL;                  % Wave number
beta = (k*b/(2*D)) .* x;      % Beta parameter

% CALCULATIONS --------------------------------------------------
IRR = (sin(beta+eps)./(beta+eps)).^2;

IRRs = (IRR.^0.2)';     % scaled irradiance for simulated screen pattern
xs = x .* 1e3;               % screen position in mm 
xs1 = x1 * 1e3; xs2 = x2 * 1e3;

% GRAPHICS ------------------------------------------------------
figure(1)
set(gcf,'Name','Single Slit Diffraction');
set(gcf,'NumberTitle','off');
set(gcf,'PaperType','A4');
set(gcf,'Color',[1 1 1]);
graphColor = ColorCode(wL);
fs = 12;
pos = [0.1 0.1 0.33 0.48];
set(gcf,'Units','normalized');
set(gcf,'Position',pos);


% Irradiance - line graph
axes('position',[0.1 0.35 0.8 0.5]);
hplot = plot(xs,IRR);
set(hplot,'Color',graphColor,'LineWidth',2);
%plot(xs,IRR,'LineWidth',2,[1 0 0]);
axis([xs1 xs2 0 1.1])
xlabel('screen position  [mm]','fontSize',fs);
ylabel('Irradiance  [ W.m^{-2}]','fontSize',fs);

mainTitle  = sprintf('Single Slit');
pT1 = '\\lambda = ';
pT2 = num2str(wL,'%3.1e');
pT3 = ' m    b = ';
pT4 = num2str(b,'%3.1e');
pT5 = ' m ';
pT6 = '    D = ';
pT7 = num2str(D,'%3.1f');
pT8 = ' m';

pT = [pT1 pT2 pT3 pT4 pT5 pT6 pT7  pT8];
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







