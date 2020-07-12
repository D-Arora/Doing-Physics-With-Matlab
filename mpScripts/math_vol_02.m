% math_vol_01.m
% 12 may 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
clear all
close all
clc

% INPUTS =================================================================
% number of partitions
    N = 2001;     % must be odd
% xA   xB  lower and upper bounds of region to be rotated
   xA = 2;
   xB = 5;

% x1   x2  limits for plotting function y = f(x)
   x1 = 0;
   x2 = 6;

% x values for region   /   xf values for function 
   x  = linspace(xA,xB,N);
   xF = linspace(x1,x2,N);
  
% Region y   /   Function  yF
  % y  = 2 .* sqrt(x);
   %yF = 2 .* sqrt(xF);
  k = 5;
   y  = 2+(x.^4 - k .* x.^2)/4;
   yF = 2+(xF.^4 - k .* xF.^2)/4;
   %y = k .* sin(2*pi*x/2);
  %y = -0.5 .* x + 2;
  %yF = -0.5 .* xF + 2 ;
  %y = ones(1,N) .* 2;
  %yF = y;
   
   yA = y(1); yB = y(end);
   
   
% Volume calculation by disk method
fn = y.^2; a = xA; b = xB;
vol_pie = simpson1d(fn,a,b);

disp('volume/pie');
disp(vol_pie);


% Setup for animated gif
ag_name = 'ag_volume.gif';   % file name for animated gif
delay = 1;                   % A scalar value 0 to 655 (inclusive), which 
                             % specifies the delay in seconds before
                             % displaying the next image.
   

% GRAPHICS ===============================================================

figure(1)    % Region   /    Function   ---------------------------------
    col = 'b';
    fs = 9;
    set(gcf,'units','normalized','position',[0.1 0.7 0.22 0.22]);
    tx = 'x_{old}';
    ty = 'y_{old} ';
    
    xP = xF; yP = yF;
    plot(xP,yP,col,'lineWidth',2);
    xlabel(tx); ylabel(ty);
    hold on
   
    xP = x; yP = y;
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[0.8 0.8 1]);
    
    xP = [x1 x2]; yP = [2 2];
    col = 'r';
    plot(xP,yP,col,'lineWidth',2);
    text(1,2.2,'axis of rotation');
   %set(gca,'Ylim',[0 1.1 * max(y)]);
%set(gca,'Ylim',Yrange);
%set(gca,'Xtick',xValues);
    grid on;
    box on;

figure(2)   % [3D] plot -------------------------------------------------
set(gcf,'units','normalized','position',[0.35 0.7 0.22 0.22]);
[X,Y,Z] = cylinder(y,100);
z = xA + Z.*(xB-xA);
surf(z,Y,X);
%axis square
xlabel('X'); ylabel('Y'); zlabel('Z') 
shading interp
box on
view(9,12);

figure(3)   % [3D] plot --------------------------------------------------
set(gcf,'units','normalized','position',[0.6 0.7 0.22 0.22]);
[X,Y,Z] = cylinder(y,100);
z = xA + Z.*(xB-xA);
surfl(z,Y,X);
axis square
xlabel('X'); ylabel('Y'); zlabel('Z') 
shading interp
colormap(copper)
hold on
plot3([xA xB],[0 0], [0 0]);
box on
view(23,16);

figure(4)   % animation of the rotation ---------------------------------
set(gcf,'units','normalized','position',[0.1 0.4 0.22 0.22]);

% X-axis
plot3([xA xB],[0 0], [0 0],'k','linewidth',2);
%set(gca,'Xlim',[xA xB]);
%set(gca,'Ylim',[xA xB]);
%set(gca,'Zlim',[xA xB]);
xlabel('x'); ylabel('y'); zlabel('z') 
hold on
grid on
box on
view(24,8);

% Function
xP = x; yP = y; zP = zeros(1,N);
plot3(xP,yP,zP,'b','linewidth',2);

t = 0 : pi/100 : 2*pi;
xP = x(end) .* ones(1,length(t)); yP = y(end) .* cos(t); zP = y(end).*sin(t);
h_plot3 = plot3(xP,yP,zP,'r','linewidth',1);
xP = x(1) .* ones(1,length(t)); yP = y(1) .* cos(t); zP = y(1).*sin(t);
h_plot3 = plot3(xP,yP,zP,'r','linewidth',1);

% Setup for saving images (im) 
Nt = 6;
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');  %RGB to indexed images
im(1,1,1,Nt) = 0;

c = 1;
for t = 0:pi/Nt:2*pi;
xP = x; yP = y .* cos(t); zP = y.*sin(t);
h_plot3 = plot3(xP,yP,zP,'r','linewidth',1);
%plot3(xB,yP(end),zP(end),'ro');
%pause(1)
%set(gca,'Xlim',[xA xB]);
%set(gca,'Ylim',[xA xB]);
%set(gca,'Zlim',[xA xB]);
 
 c = c+1;
 f = getframe(gcf);
 im(:,:,1,c) = rgb2ind(f.cdata,map,'nodither');
end

%%
figure(11)    % Region   & Reflection   ---------------------------------
    col = 'b';
    fs = 9;
    set(gcf,'units','normalized','position',[0.5 0.4 0.22 0.22]);
    tx = 'x_{new}';
    ty = 'y_{new} ';
    
    xP = x; yP = y;
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[0.8 0.8 1]);
    hold on
    xP = xF; yP = yF;
    plot(xP,yP,col,'lineWidth',2);
    xlabel(tx); ylabel(ty);
    hold on
   
      
% Reflection
      col = [1 0.8 0.8];
      xP = x; yP = -y;
      h_area = area(xP,yP);
      set(h_area,'FaceColor',col);
      
      col = 'r';
      xP = xF; yP = -yF;
      plot(xP,yP,col,'lineWidth',1);
      xlabel(tx); ylabel(ty);
         
% LIMITS
    xP = [xA xA]; yP = [0 yA];
    col = 'k';
    plot(xP,yP,col,'lineWidth',1);  
    xP = [xB xB]; yP = [0 yB];
    
    plot(xP,yP,col,'lineWidth',1);  
    
    xP = [0 xA]; yP = [yA yA];
    plot(xP,yP,col,'lineWidth',1);  
    
    xP = [0 xB]; yP = [yB yB];
    plot(xP,yP,col,'lineWidth',1);  
    
        
   %set(gca,'Ylim',[0 1.1 * max(y)]);
%set(gca,'Ylim',Yrange);
%set(gca,'Xtick',xValues);
title('Volume V_2');
    grid on;
    box on;

%%

%  SAVE ANIMATED GIF ======================================================
% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously

%imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);




