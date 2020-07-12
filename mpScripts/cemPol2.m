% cemPol2.m

% Documentation 
% http://www.physics.usyd.edu.au/teach_res/mp/doc/cemPol1.htm

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

% DOING PHYSICS WITH MATLAB 
%    http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm



% 3D plot showing two sinusoidally varying electric fields, 45 degrees out 
%   of phase.  Show resultant field from their superposition.

% waves E = E1cos(kz-wt)x^ + E2cos(kz-wt+phi)y^
% Define colours

clear all; close all; clc;

colour_darkgreen = [84 150 81] ./ 255;
%don't use normal green for projector!

vidObj = VideoWriter('avi_3DvecPolarization2.avi');
vidObj.Quality = 60;
vidObj.FrameRate = 7;
open(vidObj);

az = 110; % initial azimuth and elevation for 3D plot
el = 12;
d_az = 0; d_el = 0; % change in az and el for every plot (i.e. non zero
% values cause camera angle to change through plot

% tdisp = 0; %time displacement vector, to show progression of vectors with time.

% time
dwt = 0.2;
wt = 0:dwt:3*pi; % vetor de tempo
nt = numel(wt);

%  space
dkz = 0.01;
kz = 0:dkz:3*pi;
nz = numel(kz);

% initial phase
phi = +pi/4; %phase angle between cosines.


h = figure; set(gcf, 'Color','white');
xlen = 800; ylen=700;
h.Position=[0,0,xlen,ylen];

sx = zeros(nt,nz);
sy = zeros(nt,nz);
Ex = 1; Ey = 1;

% create electric fields in x and y planes
for t = 1:nt
    sx(t,:) = Ex*cos(kz - wt(t));
    sy(t,:) = Ey*cos(kz+phi- wt(t));
end

for t = 1:nt
    clf(h); %clear figure
    
    fprintf('frame %d of %d\n',t,nt);
    
    grid on
    daspect([1 1 1])
    set(gca,'XTickLabel',' ')
    set(gca,'YTickLabel',' ')
    set(gca,'ZTickLabel',' ')
    box on
    
    az = az + d_az; %increment az and el (view angle of 3D plot) if d_az d_el not 0.
    el = el + d_el;
    
    view(az, el);
    for j = 1:20:numel(kz)
        %axis(axis) %supresses warning with ARROW.
        axis([-1 ceil(max(kz)) -1 1 -1  1]);
        arrow([kz(j),0,0],[kz(j),0,sx(t,j)], 'EdgeColor','b','FaceColor','b');
        arrow([kz(j),0,0],[kz(j),sy(t,j),0], 'EdgeColor','r','FaceColor','r');
        arrow([kz(j),0,0],[kz(j),sy(t,j),sx(t,j)], 'EdgeColor',colour_darkgreen,'FaceColor',colour_darkgreen);
    end
    text(-2,-5,0.5,sprintf('\\bfE\\rm(0,0,z,t) = \\bfi\\rm cos(k z - \\omega t) + \\bfj\\rm cos(k z - \\omega t  +  \\phi) \n \n    \\phi = %.2f \\pi ',phi/pi),'fontsize',14);
    writeVideo(vidObj, getframe(gcf,[0 0 ,xlen,ylen]));
    
end
close(h); % close figure and video objects
close(vidObj);

