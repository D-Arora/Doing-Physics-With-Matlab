function integral = simpson2d(f,ax,bx,ay,by)

%num must be odd
%1 4 2 4 ...2 4 1

num = length(f);
hx = (bx-ax)/(num-1); hy = (by-ay)/(num-1);
h = hx * hy / 9;

% evaluates two dimension Simpson coefficients ---------------------------
sc = 2*ones(num,1);
sc(2:2:num-1) = 4;
sc(1) = 1;
sc(num) = 1;

scx = meshgrid(sc,sc);
scxy = ones(num,num);
scxy(2:2:num-1,:) = scx(2:2:num-1,:)*sc(2);
scxy(3:2:num-2,:) = scx(3:2:num-2,:)*sc(3);
scxy(1,:) = sc';
scxy(num,:) = sc';

% evaluates integral -----------------------------------------------------
integral = h * sum(sum(scxy .* f));

% simpson2d.m ************************************************************
% 250814

% Numerical integration: (2D)  - surface integrals

% Users functiond
%     simpsonxy_coeff.m

% 23 aug 2014
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/optics/optics_home.htm

% ************************************************************************