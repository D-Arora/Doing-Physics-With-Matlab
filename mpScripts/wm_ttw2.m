% wm_ttw2.m

% [3D] view of the wave function for a travelling wave 
%    propagating in the Z direction

% John A Sims 
% email: john.sims@ufabc.edu.br
% Universidade Federal do ABC 
% http://www.ufabc.edu.br/
% Biomedical Engineering Department,
% Centro de Engenharia, Modelagem e Ciencias Sociais Aplicadas.
% (Centre for Engineering, modelling and applied social sciences)
% Rua Arcturus (Jd Antares)
% Anchieta
% 09606070 - Sao Bernardo do Campo, SP - Brasil

% Ian Cooper
% email: ian.cooper@sydney.edu.au
% School of Physics, University of Sydney
% 180127
% DOING PHYSICS WITH MATLAB 
%    https://d-arora.github.io/Doing-Physics-With-Matlab/
%    http://www.physics.usyd.edu.au/teach_res/mp/doc/wm_ttw.htm


clear all
close all
clc


% =======================================================================
% TRAVELLING WAVE INPUT PARAMETERS
%   S.I. units for all quantities
%   Default values given in brakets
% =======================================================================
   A = 10;       %  amplitude  (10)
   L = 25;       %  wavelength  (25)
   T = 20;       %  period  (20)
   nT = 3;       %  no. periods for animation  (4)
   nX = 4;       %  no. wavelengths displayed (3)
   Nx = 1000;    %  Number of grid X grid points (2076)
   Nt = 500;     %  Number of time steps (200)
   phi = 0;         % phase angle [0 rad]

   
% =======================================================================
%   CALCULATIONS
% =======================================================================
   w = 2*pi / T;    %  angular frequency 
   k = 2*pi / L;    %  wave number k
   tMax = nT * T;   %  max time interval for animation
   zMax = nX * L;   %  max Z distance
  
   t = linspace(0,tMax,Nt);   % grid points for time 
   z = linspace(0,zMax,Nx);   % grid points for Z axis
   
    
%  Wave function s(z,t) ---------------------------------------------
%    Calculate all data and store it in S, rows are z, cols t. 
%    Later, phase relationship between all particles in an instant 
%    of time (i.e spatial domain information) will be taken from , 
%    rows of S and time (temporal domain) for a single particle 
%    taken from columns.
   s = zeros(Nx,Nt); 
   
% Time loop   
   for c = 1 : Nt
     s(:,c) = A .* sin(k.*z - w*t(c) + phi);
   end

% Wavefunction at a constant z and a constant t   
     Pt = find(t > 20,1)-1;
     Pz = find(z > 60,1)-1;
     st = s(Pt,:);
     
     sz = s(:,Pt);
   
% =======================================================================     
%   GRAPHICS 
% =======================================================================
figure(1)
pos = [0.05 0.1 0.3 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   [zz, tt] = meshgrid(z,t);
   mesh(zz,tt,s');
   colorbar;
   xlabel('z');
   ylabel('t');
   zlabel('s');
   view(-13,75);
   set(gca,'fontsize',12);
   tm1 = 's  =  A sin(2\pi (z / \lambda - t / T) + \phi)   \phi  =  ';
   tm2 = num2str(phi/pi, '%3.1f');
   tm3 = ' \pi rad';
   tm = [tm1 tm2 tm3];
   title(tm)
   
figure(2)
pos = [0.37 0.1 0.3 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   xP = z; yP = sz;
   plot(xP,yP,'linewidth',2);
  
   xlabel('z');
   ylabel('s');
   tm1 = 'wavefunction s at time t =  ';
   tm2 = num2str(t(Pt), '%3.1f \n');
   tm3 = '    \lambda  =  ';
   tm4 = num2str(L,'%2.1f \n');
   tm = [tm1 tm2 tm3 tm4];
   title(tm);
   grid on
   set(gca,'xTick',0:10:100);
   set(gca,'fontsize',12);
   
 figure(3)
   pos = [0.69 0.1 0.3 0.3];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w'); 
   xP = t; yP = st;
   plot(xP,yP,'linewidth',2);
   xlabel('t');
   ylabel('s');
   tm1 = 'wavefunction s at position  z  =  ';
   tm2 = num2str(z(Pz), '%3.1f \n');
   tm3 = '    T  =  ';
   tm4 = num2str(T,'%2.1f \n');
   tm = [tm1 tm2 tm3 tm4];
   title(tm);
   grid on
   set(gca,'fontsize',12);