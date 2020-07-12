% wav_shm_sine.m

clear all
clc
close all


xF = 20;
yF = 50;
xL = 1200;
yL = 800;

% Input Initial Data
ym = 10; T = 5; phi = 0; tMax = 30; Nt = 500;
boxA = ym; boxB = T; boxC = phi; boxD = tMax;

t = linspace(0,tMax,Nt);
y = ym .* sin(2*pi*t/T + phi);

f = 1/T;  w = 2*pi*f;

% Main figure window -----------------------------------------------------
f1 = figure('Color',[0.95 0.9 0.9],'Name','SHM', ...
     'NumberTitle','off','Position',[xF yF xL yL]);

 % heading ---------------------------------------------------------------
 pos = [250 750 600 30];
 colorBG = [0.95 0.9 0.9];
 colorFG = [0 0.5 1];
 fs = 18;
 textD = 'Simple Harmonic Motion - Sine Function';
 t1 = uicontrol(gcf,'Style','text','Position',pos, ...
      'String',textD,'FontSize',fs, ...
      'HorizontalAlignment','center','FontWeight','bold', ...
      'BackgroundColor',colorBG,'ForegroundColor',colorFG);

  
% INPUT WINDOW+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 plot1 = subplot('Position',[0.01 0.2 0.2 0.7]);
 set(gca,'Xlim',[0 10]);
 set(gca,'Ylim',[0 10]);
%  %set(gca,'OuterPosition',[0.26 0.6 0.7 0.3]);
text(0,9,'amplitude A','FontSize',12');
text(0,8.5,'0 to 10  [m]','FontSize',12');
text(0,7,'period  T  [s]','FontSize',12');
text(0,5,'initial phase angle','FontSize',12');
text(0,4.5,'      \phi ','FontSize',12');
text(0,4,'0 to 2 \pi  [rad]','FontSize',12');
text(0,2,'max display time  [s]','FontSize',12');
%  tm1 = 'frequency  f = ';
%  tm2 = num2str(f,'%3.3f\n');
%  tm3 = '  Hz';
%  tm = [tm1 tm2 tm3];
%  text(2,3,tm,'FontSize',12');
%  
%  tm1 = 'angular frequency  \omega = '; tm2 = num2str(w,'%3.3f');
%  tm3 = '  rad/s';
%  tm = [tm1 tm2 tm3];
%  text(2,1,tm,'FontSize',12');
axis off


% box 1 A ----------------------------------------------------------------- 
 pos = [220 630 100 30];
 colorBG = [1 1 1];
 colorFG = [0 0 0];
 fs = 14;
 Edit_A = uicontrol(gcf,'Style','edit','Position',pos, ...
          'String',boxA,'FontSize',fs,'BackgroundColor',colorBG, ...
          'Callback','boxA = str2num(get(Edit_A,''String''));'); 

  % box 2 B -----------------------------------------------------------------
  pos = [220 540 100 30];
  colorBG = [1 1 1];
  colorFG = [0 0 0];
  fs = 14;
  Edit_B = uicontrol(gcf,'Style','edit','Position',pos, ...
          'String',boxB,'FontSize',fs,'BackgroundColor',colorBG, ...
          'Callback','boxB = str2num(get(Edit_B,''String''));');
  
% box 3 C -----------------------------------------------------------------
  pos = [220 390 100 30];
  colorBG = [1 1 1];
  colorFG = [0 0 0];
  fs = 14;
  Edit_C = uicontrol(gcf,'Style','edit','Position',pos, ...
          'String',boxC,'FontSize',fs,'BackgroundColor',colorBG, ...
          'Callback','boxC = str2num(get(Edit_C,''String''));');
  
% box 4 D -----------------------------------------------------------------
  pos = [220 250 100 30];
  colorBG = [1 1 1];
  colorFG = [0 0 0];
  fs = 14;
  Edit_D = uicontrol(gcf,'Style','edit','Position',pos, ...
          'String',boxD,'FontSize',fs,'BackgroundColor',colorBG, ...
          'Callback','boxD = str2num(get(Edit_D,''String''));');
       
  
 % PUSHBUTTONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
pushbutton_run = uicontrol(gcf,'Style','pushbutton','Position',...
     [50 50 100 40], 'FontSize',12,'FontWeight','bold', 'String','RUN', ...
     'CallBack', 'wav_shm_sine_cal');
 
 
 pushbutton_close = uicontrol(gcf,'Style','pushbutton','Position',...
     [200 50 100 40], 'FontSize',12,'FontWeight','bold', 'String','CLOSE', ...
     'CallBack', 'close');
  
 % PLOTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 plot1 = subplot('Position',[0.4 0.1 0.5 0.4]);
    xP = t; yP = y;
    plot(xP,yP,'k','lineWidth',2);
    axis on; grid on;
    xlabel('time  t  [s]','FontSize',12');
    ylabel('y  [m]','FontSize',12');
    set(gca,'Ylim',[-10 10]);
 % -----------------------------------------------------------------------
    
 plot1 = subplot('Position',[0.4 0.6 0.5 0.3]);
 set(gca,'Xlim',[0 10]);
 set(gca,'Ylim',[0 10]);
 %set(gca,'OuterPosition',[0.26 0.6 0.7 0.3]);
 text(2,9,'y = A sin(2\pi t / T + \phi)','FontSize',12');
 text(2,7,'y = A sin(2\pi f t + \phi)','FontSize',12');
 text(2,5,'y = A sin(\omega t + \phi)','FontSize',12');
 
 tm1 = 'frequency  f = ';
 tm2 = num2str(f,'%3.3f\n');
 tm3 = '  Hz';
 tm = [tm1 tm2 tm3];
 text(2,3,tm,'FontSize',12');
 
 tm1 = 'angular frequency  \omega = '; tm2 = num2str(w,'%3.3f');
 tm3 = '  rad/s';
 tm = [tm1 tm2 tm3];
 text(2,1,tm,'FontSize',12');
 axis off
 

     