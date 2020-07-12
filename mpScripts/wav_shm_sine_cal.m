% wav_shm_sine_cal.m
% CallBack mscript for wav_shm_sine.m

% Reads values from input boxes -----------------------------------------
ym = boxA;
T  = boxB;
phi = boxC;
tMax = boxD;

t = linspace(0,tMax,Nt);
y = ym .* sin(2*pi*t/T + phi);

f = 1 / T; w = 2*pi*f;

       
 % PLOT ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 plot1 = subplot('Position',[0.4 0.1 0.5 0.4]);
    xP = t; yP = y;
    plot(xP,yP,'k','lineWidth',2);
    axis on; grid on;
    xlabel('time  t  [s]','FontSize',12');
    ylabel('y  [m]','FontSize',12');
    set(gca,'Ylim',[-10 10]);

 % Output parameters -----------------------------------------------------
    
 plot1 = subplot('Position',[0.4 0.6 0.5 0.3]);
 set(gca,'Xlim',[0 10]);
 set(gca,'Ylim',[0 10]);
 %set(gca,'OuterPosition',[0.26 0.6 0.7 0.3]);
 text(2,9,'y = A sin(2\pi t / T + \phi)','FontSize',12');
 text(2,7,'y = A sin(2\pi f t + \phi)','FontSize',12');
 text(2,5,'y = A sin(\omega t + \phi)','FontSize',12');
 
 colorBG = [0.95 0.9 0.9];
 tmA = 'yyyyyyyyyyyyyyyyyyyyyyyyyyyyyy';
 text_h = text(2,3,tmA,'FontSize',16','color',colorBG,'EdgeColor',colorBG, ...
     'BackgroundColor',colorBG);
 
 tm1 = 'frequency  f = ';
 tm2 = num2str(f,'%3.3f\n');
 tm3 = '  Hz';
 tm = [tm1 tm2 tm3];
 text(2,3,tm,'FontSize',12');
 
 text_h = text(2,1,tmA,'FontSize',16','color',colorBG,'EdgeColor',colorBG, ...
     'BackgroundColor',colorBG);
 tm1 = 'angular frequency  \omega = '; tm2 = num2str(w,'%3.3f');
 tm3 = '  rad/s';
 tm = [tm1 tm2 tm3];
 text(2,1,tm,'FontSize',12');
 axis off
     
       
       
