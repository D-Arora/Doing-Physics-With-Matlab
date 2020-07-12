% cg001.m


clear 
close all
clc

% =======================================================================
% INPUT: Force Field --> acceleration
% flagF = 1   a = (0,0)  
% flagF = 2   a = (0,-10)    gravitational force
% flagF = 3   a = (50,0)/m   electric force
  flagF = 1;

% =======================================================================   
% Setup for saving images animated gif file and AVI files
% =======================================================================
   flagA = 0;
   
%  ANIMATED GIF:
     ag_name = 'ag_cg003.gif';
%  Delay in seconds before displaying the next image  
    delay = 0.02;  
%  Frame counter start
    nt = 1;

%  AVI: open writer object for avi.
   aviName = 'am_cg003.avi';
   aviObj = VideoWriter(aviName);
   aviObj.FrameRate = 10;  % from 30
   aviObj.Quality = 70;    % from 75
   open(aviObj);

% ======================================================================   
L = 100;
tMax = 5;
tSteps = 500;
t = linspace(0,tMax,tSteps);
dt = t(2) - t(1);
m = 10:5:30;

z1 = zeros(tSteps,1) ;
v1 = zeros(tSteps,1) ; 
z1(1) = 75 + 1i*50;
v1(1) = -10 + 1i*10;

z2 = zeros(tSteps,1) ;
v2 = zeros(tSteps,1) ; 
z2(1) = -75 + 1i*50;
v2(1) = +10 + 1i*10;

z3 = zeros(tSteps,1) ;
v3 = zeros(tSteps,1) ; 
z3(1) = -25 + 1i*0;
v3(1) = -10 + 1i*20;

z4 = zeros(tSteps,1) ;
v4 = zeros(tSteps,1) ; 
z4(1) = 50 - 1i*0;
v4(1) = 20 - 1i*30;

z5 = zeros(tSteps,1) ;
v5 = zeros(tSteps,1) ; 
z5(1) = 0 - 1i*0;
v5(1) = 0 + 1i*0;



for c = 2 : tSteps
[z1(c), v1(c)]  = posVel(z1(c-1),v1(c-1),dt,m(1),flagF); 
[z2(c), v2(c)]  = posVel(z2(c-1),v2(c-1),dt,m(2),flagF);
[z3(c), v3(c)]  = posVel(z3(c-1),v3(c-1),dt,m(3),flagF);
[z4(c), v4(c)]  = posVel(z4(c-1),v4(c-1),dt,m(4),flagF);
[z5(c), v5(c)]  = posVel(z5(c-1),v5(c-1),dt,m(5),flagF);
end

%%
figure(1)   % 11111111111111111111111111111111111111111111111111111111111
   set(gcf,'units','normalized','position',[0.05 0.05 0.33 0.42]);
   tm1 = 't  =  ';
   tm3 = '  s';
   
   for c = 1 : tSteps
       tm2 = num2str(t(c),'%3.2f \n');
       tm = [tm1 tm2 tm3];
       xP = real(z1(c)); yP = imag(z1(c));
       hPlot = plot(xP,yP,'o');
       set(hPlot,'MarkerSize',10,'MarkerFaceColor',[1 0 0 ], 'MarkerEdgecolor',[1 0 0])
       hold on
       xP = real(z2(c)); yP = imag(z2(c));
       hPlot = plot(xP,yP,'o');
       set(hPlot,'MarkerSize',15,'MarkerFaceColor',[0 0 1 ], 'MarkerEdgecolor',[0 0 1])
       xP = real(z3(c)); yP = imag(z3(c));
       hPlot = plot(xP,yP,'o');
       set(hPlot,'MarkerSize',20,'MarkerFaceColor',[1 0 1 ], 'MarkerEdgecolor',[1 0 1])
       xP = real(z4(c)); yP = imag(z4(c));
       hPlot = plot(xP,yP,'o');
       set(hPlot,'MarkerSize',25,'MarkerFaceColor',[0 0 0], 'MarkerEdgecolor',[0 0 0])
       xP = real(z5(c)); yP = imag(z5(c));
       hPlot = plot(xP,yP,'o');
       set(hPlot,'MarkerSize',30,'MarkerFaceColor',[1 1 0], 'MarkerEdgecolor',[0 0 0])
       
       set(gca,'xLim',[-L L])
       set(gca,'yLim',[-L L])
       axis square
       grid on
       xlabel('x  [ m ]')
       ylabel('y  [ m ]')
       title(tm)
       set(gca,'fontsize',14)
       pause(0.01) 
       hold off
       
       if flagA == 1
         frame = getframe(1);
         im = frame2im(frame);
         [imind,cm] = rgb2ind(im,256);
       %  On the first loop, create the file. In subsequent loops, append.
          if nt == 1
           imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'loopcount',inf);
          else
          imwrite(imind,cm,ag_name,'gif','DelayTime',delay,'writemode','append');
          end
        nt = nt+1;
          writeVideo(aviObj,frame);
       end
       
   end
        if flagA == 1
         close(aviObj);
        end
   
function [z, v] = posVel(z1,v1,dt,m, flagF)
  
 L = 100;
 if flagF == 1; a = 0;              end      % zero acceleration
 if flagF == 2; a = 0 - 1i*10;      end      % gravitational field
 if flagF == 3; a = (50 - 1i*0)/m;  end      % horizontal force -->
 
  
v = v1 + a*dt;
  if real(z1) > L;  v = -real(v) + 1i*imag(v); end
  if real(z1) < -L; v = -real(v) + 1i*imag(v); end
  if imag(z1) > L;  v = real(v)  - 1i*imag(v); end
  if imag(z1) < -L; v = real(v)  - 1i*imag(v); end
  z = z1 + v*dt + 0.5 * a * dt^2;

end
 


