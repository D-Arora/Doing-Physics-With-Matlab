%empulse02.m
%b01/mat/em
%11 sep 01

%FDTD simulation of a pulse in free space
%A Gaussian pulse originates in the centre and travels outwards in z direction
%electric field Ex
%magnetic field Hy
%absorbing boundary conditions

tic
close all
clear all
clc

ag_name = 'ag_empulse01a.gif';   % file name for animated gif
delay = 0.1;                   % A scalar value 0 to 655 (inclusive), which 
                             % specifies the delay in seconds before
                             % displaying the next imag

% Define z axis   --------------------------------------------------
zmin = 0;
zmax = 200;
nz = 400;
dz = (zmax - zmin)/(nz-1);
z = zmin : dz : zmax;

% Define time t    ------------------------------------------------
nt = 220;    %number of time steps

% Initialise arrays    ----------------------------------------------
E = zeros(nz,1);
H = zeros(nz,1);

% Gaussian pulse   ---------------------------------------------------
%zc = zmax/2;   %z position of pulse centre 
zc = nz/2;
s = 12;       %width of pulse
A = 1;        %pulse height
to = 40;

% Graphics ------------------------------------------------------------
figure(1)
subplot(2,1,1);
plot(z,E,'LineWidth',2);
xlabel('position  z')
ylabel('electric field   Ex');
set(gca,'Ylim',[-1.2 1.2]);
grid on

subplot(2,1,2);
plot(z,H);
xlabel('position  z')    
ylabel('magnetic field   Hy');
grid on

% Setup for saving images (im) 
Nt = nt/10;
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');  %RGB to indexed images
im(1,1,1,Nt) = 0;

% Gaussian pulse
Ep = zeros(1,nt); step = zeros(1,nt);
for ct = 1 : nt
Ep(ct) = A*exp(-0.5*((to-ct)/s)^2);
step(ct) = ct;

end
figure(2);
plot(step,Ep,'linewidth',2);
xlabel('time step');
ylabel('E_x(at source point)');
title('Hard Source');


% Calculations   ---------------------------------------------------------
c = 1;
for ct = 1 : nt
   for cz = 2 : nz-1
      E(cz) = E(cz) + 0.5*(H(cz-1)-H(cz));
      pulse = A*exp(-0.5*((to-ct)/s)^2);
      if pulse > 0.02, E(zc) = pulse; end;
   end   % for cz
   
   for cz = 2 : nz-1
      H(cz) = H(cz) + 0.5*(E(cz) - E(cz+1));
   end   % for cz
   
   
  if mod(ct,10) == 0;
  figure(1)
  subplot(2,1,1)
  plot(z,E,'color',[0 0 1],'LineWidth',2);
  set(gca,'Ylim',[-1.2 1.2]);
  xlabel('position  z (m)')
  ylabel('electric field   Ex');
  titleString = ['time step   ',num2str(ct)]; 
  title(titleString);
  
  subplot(2,1,2);
  plot(z,H,'color',[1 0 0],'lineWidth',2);
  set(gca,'Ylim',[-1.2 1.2]);
  xlabel('position  z (m)')    
  ylabel('magnetic field   Hy');

  pause(0.1);
  drawnow;
  f = getframe(gcf);
  im(:,:,1,c) = rgb2ind(f.cdata,map,'nodither');
  c = c+1;
  end   %if

end   % for ct   

%  SAVE ANIMATED GIF ======================================================
% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously

%imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);

toc


