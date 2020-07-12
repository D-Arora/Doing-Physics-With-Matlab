% math_vol_01.m
% 12 may 2015
% Ian Cooper   School of Physics   University of Sydney
% cooper@physics.usyd.edu.au
% https://d-arora.github.io/Doing-Physics-With-Matlab/
clear all
close all
clc

% INPUTS =================================================================
% number of partitions
    N = 2001;     % must be odd
% xA   xB  lower and upper bounds of region to be rotated
   xA = 2;
   xB = 4;

% x1   x2  limits for plotting function y = f(x)
   x1 = 0;
   x2 = 6;

% axis of rotation
   yR = +2;
% x values for region   /   xf values for function 
   x  = linspace(xA,xB,N);
   xF = linspace(x1,x2,N);
  
% Region y   /   Function  yF
  % y  = 2 .* sqrt(x);
   %yF = 2 .* sqrt(xF);
   %k = 5;
   %y  = 2+(x.^4 - k .* x.^2)/4;
   %yF = 2+(xF.^4 - k .* xF.^2)/4;
   %y = k .* sin(2*pi*x/2);
  y = 0.5 .* x;
  yF = 0.5 .* xF;
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
    tx = 'x';
    ty = 'y';
    
    hold on
    
    xP = x; yP = y;
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[0.8 0.8 1]);
    
    xP = xF; yP = yF;
    plot(xP,yP,col,'lineWidth',2);
    xlabel(tx); ylabel(ty);
    title('Function y = f(x) and Region');     
% axis of rotation
    xP = [x1 x2]; yP = [yR yR];
    col = 'r';
    plot(xP,yP,col,'lineWidth',2);
    text(0.5,2.4,'axis of rotation','Color','r');
    
    
   set(gca,'Ylim',[-3 3]);
   %set(gca,'Ylim',Yrange);
   %set(gca,'Xtick',xValues);
   grid on;
   box on;

figure(2)    % Region   /    Function   ---------------------------------
    col = 'b';
    fs = 9;
    set(gcf,'units','normalized','position',[0.35 0.7 0.22 0.22]);
    tx = 'x';
    ty = 'y';
        hold on
        

        
% region - reflection
    xP = [xA xB]; yP = [2*yR 2*yR];
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[0.8 0.8 1]);
    
    xP = x; yP = -y+2*yR;
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[1 1 1]); 
           
    xP = x; yP = y;
    h_area = area(xP,yP);
    set(h_area,'FaceColor',[0.8 0.8 1]);
    
% function
    xP = xF; yP = yF;
    plot(xP,yP,col,'lineWidth',2);    

% function - reflection
    col = 'm';
    xP = xF; yP = -yF+2*yR;
    plot(xP,yP,col,'lineWidth',2);
    xlabel(tx);
    h = ylabel(ty);
    hold on
    
% axis of rotation
    xP = [x1 x2]; yP = [yR yR];
    col = 'r';
    plot(xP,yP,col,'lineWidth',2);
    text(0.1,2.5,'axis of rotation','color','r');

% y = yR;
    col = 'r'; LW = 2;
    xP = [2 4]; yP = [2*yR 2*yR];
    plot(xP,yP,col,'lineWidth',LW);
    
    
    xlabel(tx);
    ylabel(ty,'Rotation',0); 
    title('Rotated Region in XY plane');
    
   set(gca,'Ylim',[-2  4.2]);
   set(gca,'Xlim',[0 4.2]);
   %set(gca,'Xtick',xValues);
   %grid on;
   box on;    
    
figure(6)   % [3D] plot -------------------------------------------------
set(gcf,'units','normalized','position',[0.35 0.4 0.22 0.22]);
[X,Y,Z] = cylinder(y-yR,100);
z = xA + Z.*(xB-xA);
surf(z,Y+yR,X);
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z') 
shading interp
box on
view(35,10);
axis off
set(gca,'Xlim',[0 4]);
set(gca,'Ylim',[-2 4]);
set(gca,'Zlim',[-2 4]);

figure(3)   % [3D] plot --------------------------------------------------
set(gcf,'units','normalized','position',[0.6 0.7 0.22 0.22]);
Rad = ones(1,N) .* yR;
[X,Y,Z] = cylinder(Rad,100);
z = xA + Z.*(xB-xA);
surf(z,Y+yR,X);
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z') 
shading interp
box on
view(35,10);
axis off
set(gca,'Xlim',[0 4]);
set(gca,'Ylim',[-2 4]);
set(gca,'Zlim',[-2 4]);


% 
figure(4)   % animation of the rotation ---------------------------------
set(gcf,'units','normalized','position',[0.1 0.4 0.22 0.22]);

% OUTTER SURFACE
% X Y Z axes
   col = 'k';  LW = 1;
   plot3([0 6], [0 0], [0 0],col,'linewidth',LW);
   hold on
   plot3([0 0],[-6 6],[0 0],col,'linewidth',LW);
   plot3([0 0], [0 0], [-6 6],col,'linewidth',LW);
   xlabel('x'); ylabel('y'); zlabel('z') 
   hold on
   grid on
   box on
   view(42,18);

% Function
   xP = x; yP = y; zP = zeros(1,N);
   plot3(xP,yP,zP,'b','linewidth',2);

% Circles radius R
   col = 'b'; LW = 1;
   R = -y + yR;

   t = 0 : pi/100 : 2*pi;
   xP = x(1) .* ones(1,length(t));
   yP = yR + R(1).* cos(t);
   zP = R(1) .* sin(t);
   plot3(xP,yP,zP,col,'linewidth',LW);

   xP = x(end) .* ones(1,length(t));
   yP = yR + R(end) .* cos(t);
   zP = R(end) .* sin(t);
   plot3(xP,yP,zP,col,'linewidth',LW);

% Setup for saving images (im) 
%Nt = 6;
%f = getframe(gcf);
%[im,map] = rgb2ind(f.cdata,256,'nodither');  %RGB to indexed images
%im(1,1,1,Nt) = 0;

% Profile
Nt = 25;
c = 1;
for t = 0:pi/Nt:2*pi;
   xP = x;
   yP = yR + R .* cos(t);
   zP = R.*sin(t);
   plot3(xP,yP,zP,col,'linewidth',LW);

%plot3(xB,yP(end),zP(end),'ro');
%pause(1)
%set(gca,'Xlim',[xA xB]);
%set(gca,'Ylim',[xA xB]);
%set(gca,'Zlim',[xA xB]);
 
 %c = c+1;
 %f = getframe(gcf);
 %im(:,:,1,c) = rgb2ind(f.cdata,map,'nodither');
end

% INNER SURFACE

% Function
   xP = x; yP = zeros(1,N); zP = zeros(1,N);
   plot3(xP,yP,zP,'m','linewidth',2);

% Circles radius R
   col = 'r'; LW = 1;
   R = yP - yR;

   t = 0 : pi/100 : 2*pi;
   xP = x(1) .* ones(1,length(t));
   yP = yR + R(1).* cos(t);
   zP = R(1) .* sin(t);
   plot3(xP,yP,zP,col,'linewidth',LW);

   xP = x(end) .* ones(1,length(t));
   yP = yR + R(end) .* cos(t);
   zP = R(end) .* sin(t);
   h_plot3 = plot3(xP,yP,zP,col,'linewidth',1);

% Setup for saving images (im) 
%Nt = 6;
%f = getframe(gcf);
%[im,map] = rgb2ind(f.cdata,256,'nodither');  %RGB to indexed images
%im(1,1,1,Nt) = 0;

% Profile
Nt = 25;
c = 1;
for t = 0:pi/Nt:2*pi;
   xP = x;
   yP = yR + R .* cos(t);
   zP = R.*sin(t);
h_plot3 = plot3(xP,yP,zP,col,'linewidth',1);

%plot3(xB,yP(end),zP(end),'ro');
%pause(1)
%set(gca,'Xlim',[xA xB]);
%set(gca,'Ylim',[-4 5]);
%set(gca,'Zlim',[xA xB]);
 
 %c = c+1;
 %f = getframe(gcf);
 %im(:,:,1,c) = rgb2ind(f.cdata,map,'nodither');
end

set(gca,'Xlim',[0 4]);
set(gca,'Ylim',[-1 4]);
set(gca,'Zlim',[-1 4]);
set(gca,'Xtick',[0 :2: 4]);
set(gca,'Ytick',[0 :2:  4]);
set(gca,'Ztick',[0 :2: 4]);
axis equal
%
% figure(11)    % Region   & Reflection   ---------------------------------
%     col = 'b';
%     fs = 9;
%     set(gcf,'units','normalized','position',[0.5 0.4 0.22 0.22]);
%     tx = 'x_{new}';
%     ty = 'y_{new} ';
%     
%     xP = x; yP = y;
%     h_area = area(xP,yP);
%     set(h_area,'FaceColor',[0.8 0.8 1]);
%     hold on
%     xP = xF; yP = yF;
%     plot(xP,yP,col,'lineWidth',2);
%     xlabel(tx); ylabel(ty);
%     hold on
%    
%       
% % Reflection
%       col = [1 0.8 0.8];
%       xP = x; yP = -y;
%       h_area = area(xP,yP);
%       set(h_area,'FaceColor',col);
%       
%       col = 'r';
%       xP = xF; yP = -yF;
%       plot(xP,yP,col,'lineWidth',1);
%       xlabel(tx); ylabel(ty);
%          
% % LIMITS
%     xP = [xA xA]; yP = [0 yA];
%     col = 'k';
%     plot(xP,yP,col,'lineWidth',1);  
%     xP = [xB xB]; yP = [0 yB];
%     
%     plot(xP,yP,col,'lineWidth',1);  
%     
%     xP = [0 xA]; yP = [yA yA];
%     plot(xP,yP,col,'lineWidth',1);  
%     
%     xP = [0 xB]; yP = [yB yB];
%     plot(xP,yP,col,'lineWidth',1);  
%     
%         
%    %set(gca,'Ylim',[0 1.1 * max(y)]);
% %set(gca,'Ylim',Yrange);
% %set(gca,'Xtick',xValues);
% title('Volume V_2');
%     grid on;
%     box on;
% 
% %%
% 
% %  SAVE ANIMATED GIF ======================================================
% % im - images to be saved
% % map - color map for images
% % ag_name - file name for animated gif
% % DelayTime - time delay in seconds between viewing images
% % LoopCount - animated gif will continuously
% 
% %imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);
% 
% 
% 
% 
