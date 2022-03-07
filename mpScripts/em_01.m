% em_01.m

% Vector addtion and subtraction of three vectors in [2D]

% DOING PHYSICS WITH MATLAB
%    https://d-arora.github.io/Doing-Physics-With-Matlab/

% DOWNLOAD Scripts from
%   https://github.com/D-Arora/Doing-Physics-With-Matlab/tree/master/mpScripts
%   https://drive.google.com/drive/u/3/folders/1j09aAhfrVYpiMavajrgSvUMc89ksF9Jb

% 220225  Matlab Version R2021b


close all
clear
clc


% INPUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Vectors V [Vx Vy Vz = 0]
 A = [3 8 0];
 B = [-3 4 0];
 C  =[ -7 -8 0];

% Arrow parameters
  colArrow = [1 0 0; 0 0 1; 1 0 1; 0 0 0];
  LW = 2;
  L = 1;
  W = 0.5;

% CALCULATIONS --------------------------------------------------------

% Resultant vector R
  R = A+B+C;

% Input parameters for arrows in Figure Window  
V = zeros(4,3);
V(1,:) = A; V(2,:) = B; V(3,:) = C; V(4,:) = R;

% Magntiudes  mag
  mag = zeros(4,1);
  for c = 1 : 4
     mag(c) = norm(V(c,:));
  end

% Direction angle   DA
  theta = atan2(V(:,2), V(:,1));

% Components x, y, z
  x = mag.*cos(theta);
  y = mag.*sin(theta);
  z = zeros(4,1);
 

% GRAPHICS ============================================================
  figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.10 0.3 0.7]);
  set(gcf,'color','w');
  FS = 14;

subplot(2,1,1)
  for c = 1:4
      zT = [0 0 0];  
      magV = mag(c);
      angleV = theta(c);
      col = colArrow(c,:);
      DrawArrow(zT,magV,angleV,L,W,LW,col)
  end

  LIM = [-15 15];
  xlim(LIM);
  ylim(LIM);
  xticks(-10:5:10)
  yticks(-10:5:10)
  axis square
  box on
  grid on
  set(gca,'fontsize',FS)
  legend('A','B','C','R','Orientation','horizontal','Location', ...
     'northoutside','box','off','fontsize',10 )
  subplot(2,1,2)
      HT = zeros(4,3);
      HT(2,:) = [x(1) y(1) 0];
      HT(3,:) = [x(1)+x(2) y(1)+y(2) 0]; 
      
  for c = 1:4
      zT = HT(c,1)+1i*HT(c,2);
      magV = mag(c);
      angleV = theta(c);
      col = colArrow(c,:);
      DrawArrow(zT,magV,angleV,L,W,LW,col)
  end

  LIM = [-15 15];
  xlim(LIM);
  ylim(LIM);
  xticks(-10:5:10)
  yticks(-10:5:10)
  axis square
  box on
  grid on
  set(gca,'fontsize',FS)
   legend('A','B','C','R','Orientation','horizontal','Location', ...
     'northoutside','box','off','fontsize',10 )

