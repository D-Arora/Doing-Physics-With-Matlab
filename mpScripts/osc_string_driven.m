% osc-string_driven.m

% Ian Cooper
% School of Physics, University of Sydney
% Doing Physics online  
% ../mphome.htm
% email:  cooper@physics.usyd.edu.au


% Vibrations of a string driven by a sinusoidal source function
% Natural frequencies of vibration and Resonance

clear all
close all
clc

disp('   ');
disp('   String driven by a constant frequency source at x = 0');
disp('   ');
disp('   Natural frequencies of vibration: f = 625, 1250, 1875, 2500, 3125, ...  Hz')
disp('    ');

fD = input('Input driving frequency [Hz]  ');



Nx = 200;                   % no. x values
Nt = 10000;                  % no. time increments
Nstep = 25;                 % animation step
Nf = 1:5;                   % number of harmonics
NT = 5;                     % no. periods for animation
L = 0.8;                    % length of guitar strinfg [m]
FT = 400;                   % string tension [N]
mu = 4e-4;                  % linear density of string  [kg/m]

v = sqrt(FT/mu);            % speed of transverse wave along string
wL1 = 2*L;                  % fundamental wavelength   [m] 
f1 =  v / wL1;              % fundamental frequency    [Hz]
T1 = 1/f1;                  % period of fundamental
f = f1 * Nf;                % harmonic frequencies
wD = 2*pi*fD;               % angular frequency of driving force
A = 0.002;                  % amplitude of driving force

x = linspace(0,L,Nx);       % position along string
dx = x(2)-x(1);             
dt = 0.2*dx / v;            % dt < dx / v
t = 0:dt:Nt*dt;             % time

K = (v * dt / dx)^2;        % constant

y = zeros(Nt,Nx);           % initilize wave function
y(:,end) = 0;               % boundary condition at end of tube

for ct = 1:Nt
   y(ct,1) = A .* sin(wD .* t(ct));     % time variation of driving force
end

% Solving the Wave equation by FDTD method --------------------------------
for ct = 2 : Nt-1
   for cx = 2 : Nx-1
     y(ct+1,cx) = 2*y(ct,cx) - y(ct-1,cx) + K * (y(ct,cx+1) - 2* y(ct,cx) + y(ct,cx-1));
   end
end

tm1 = 'driving frequency   =   ';
tm2 = num2str(fD,'%4.0f \n');
tm3 = '   Hz';
tm = [tm1 tm2 tm3];

% Graphics ---------------------------------------------------------------
figure(1)
fs = 12;
set(gcf,'Units','Normalized');
set(gcf,'Position',[0.2 0.2 0.35 0.3]);
set(gcf,'color',[1 1 1]);
set(gca,'fontsize',12);
ymax = 0.04;
plot(x,y(1,:),'linewidth',2);
axis([0 L -ymax ymax]);
grid on
xlabel('position along string  (m)','fontsize',fs);
ylabel('string displacement   (mm)','fontsize',fs);
title(tm,'fontsize',12);

M = getframe(gcf);
[im,map] = rgb2ind(M.cdata,256,'dither');  %RGB to indexed images
im(1,1,1,1) = 0;

cc = 1;
for ct = 1: Nstep: Nt
   plot(x,y(ct,:),'linewidth',2);
   axis([0 L -ymax ymax]);
   set(gca,'Ytick',[-0.04 -0.02 0 0.02 0.04]);
   grid on
   xlabel('position along string  (m)','fontsize',12);
   ylabel('string displacement   (mm)','fontsize',12);
   title(tm,'fontsize',12)
   set(gca,'fontsize',12);
   pause(eps)
   
    M = getframe(gcf) ;
    cc = cc+1;
    im(:,:,1,cc) = rgb2ind(M.cdata,map,'dither');
   
end

% SAVE or NOT SAVE animated gif +++++++++++++++++++++++++++++++++++++++++

delay = 0;
ag_name = 'ag_driven.gif';   % file name for animated gif
imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);

