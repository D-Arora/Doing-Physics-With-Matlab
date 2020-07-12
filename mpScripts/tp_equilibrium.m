% tp_equilibrium.m
% Simulation of a N gas particles released into a box from one corner 
% Box dimensions  10 x 10,  Unit box dimensions  1 x 1
% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% ../mphome.htm
% 20140718

clear all
close all
clc

% ************************************************************************
% SETUP
% ************************************************************************

% input number of particles and time step
N = input('No. of particles (default N = 1000), N  =  ');     
disp('   ');                                            
nt = input('No. of time steps (default nt = 200),  nt  =  ');  
disp('  ');

% initial position of N particles within unit box ------------------------
R = 1;                          % max step length
P = rand(N,2);                  % initial (x,y) coordinates of N particles
x = P(:,1);
y = P(:,2);


% ************************************************************************
%  CALCULATIONS
%  GRAPHICS

% ************************************************************************ 
xp = [ 0 1 1 0 0]; yp = [0 0 1 1 0];

figure(1) % Animation of the motion of the gas molecules
set(gcf,'Units','Normalized','Position',[0.2,0.2,0.2,0.3]) 
set(gca,'fontsize',12);
%set(gcf,'Units','centimeters','Position',[1,5,10,10]);
set(gcf,'NumberTitle','off','Name','Particles in Box');

color_p = [0 0 1];  size_p = 4;

% calculate new positions of each each particle and plot -----------------     
for c = 1 : nt

% plot the position of the particles -------------------------------------
h_plot = plot(x,y,'o');
     set(h_plot,'MarkerEdgeColor',color_p,'MarkerFaceColor',color_p, ...
        'MarkerSize',size_p);
hold on
h_plot = plot(x(1),y(1),'o');
     set(h_plot,'MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[1 0 1], ...
        'MarkerSize',size_p+2);
plot(xp,yp,'r','linewidth',3');    
hold off    
    
    
axis([0 10 0 10])
axis equal
set(gca,'Xlim',[0 10]);
set(gca,'Xtick',[0:1:10]);
set(gca,'Ylim',[0 10]);
set(gca,'Ytick',[0:1:10]);
% axis off
grid on

if c ==1
M = getframe(gcf);
[im,map] = rgb2ind(M.cdata,512,'nodither');  %RGB to indexed images
im(1,1,1,nt) = 0;
end

% delay in animation for position of particles ---------------------------
if c == 1; pause(2); else pause(0.1); end;    

% calculate the number of particles in start box --------------------------
box_x = ones(1,N);
box_y = ones(1,N);
box_x(x>R) = 0;
box_y(y>R) = 0;
N_box(c) = sum(box_x .* box_y);

% increment the position of each particle --------------------------------    
dP = R.*randn(N,2);
dx = dP(:,1);
dy = dP(:,2);

% updated position of each particle---------------------------------------
x = x + dx;
y = y + dy;

% particle reflected of boundaries of the box ----------------------------
for c1 = 1:N
if x(c1) < 0; x(c1) = x(c1) + abs(dx(c1));   end;
if y(c1) < 0; y(c1) = y(c1) + abs(dy(c1));   end;
if x(c1) > 10; x(c1) = x(c1) - abs(dx(c1));  end;
if y(c1) > 10; y(c1) = y(c1) - abs(dy(c1));  end;
end  

M = getframe(gcf);
im(:,:,1,c) = rgb2ind(M.cdata,map,'nodither');
 
end

% calc. mean   std sem for no. of particles in start box ---------------
ns = floor(0.6*nt);                    % starting point for calc.
N_avg = (100/N) * mean(N_box(ns:end));
N_std = (100/N) * std(N_box(ns:end));
N_sem = N_std/sqrt(nt-ns+1);

% ======================================================================== 
figure(2)  % plot for particles in start box 
set(gcf,'NumberTitle','off','Name','Start Box');
set(gcf,'Units','centimeters','Position',[12,5,20,16]);
fs = 12; 
set(gca,'FontSize',fs);

% number of particles in unit box ----------------------------------------
subplot(2,1,1)
y_p = N_box*100/N; H = plot(y_p);

lineWidth_p = 2; set(H,'lineWidth',lineWidth_p);
    
title_x = 'time steps  (a.u)';
title_y = '% particles in sstart box';
tm1 = 'No. of particles   =   ';
tm2 = num2str(N);
tm = [tm1 tm2];
     xlabel(title_x,'FontSize',fs);
     ylabel(title_y,'FontSize',fs); 
     ht = title(tm,'FontSize',fs);

hold on

% plot red line indicating start and mean, std, sem calc. -----------------
x_p = [ns ns]; y_p = [0 100]; H = plot(x_p,y_p, 'r');

txt =  'Start box statistics'; text(ns+5,90,txt,'FontSize',fs);

txt_1 = 'mean  =  ' ; 
txt_2 = num2str(N_avg,'%0.2f');
txt   = [txt_1 txt_2];
text(ns+5,80,txt,'FontSize',fs);

txt_1 = 'std  =  ' ; 
txt_2 = num2str(N_std,'%0.2f');
txt   = [txt_1 txt_2];
text(ns+5,70,txt,'FontSize',fs);

txt_1 = 'sem  =  ' ; 
txt_2 = num2str(N_sem,'%0.2f');
txt   = [txt_1 txt_2];
%text(ns+5,60,txt,'FontSize',fs);     

% magnified view of number of particles in start box ----------------------     
subplot(2,1,2)
y_p = N_box*100/N; H = plot(y_p);

lineWidth_p = 2; set(H,'lineWidth',lineWidth_p);

set(gca,'Ylim',[0 5]);
set(gca,'Ytick',[0:0.5:5]);
    
title_x = 'time steps  (a.u)';
title_y = '% particles in start box';
tm1 = 'No. of particles   =   ';
tm2 = num2str(N);
tm = [tm1 tm2];
xlabel(title_x,'FontSize',fs);
ylabel(title_y,'FontSize',fs); 
ht = title(tm,'FontSize',fs);

     
%  SAVE ANIMATED GIF ======================================================
% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously
ag_name = 'ag_aaa.gif';
delay = 0.1;
%imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);

