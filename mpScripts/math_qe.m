% math_qe.m

% Solving Quadratic Equations - GUI
% Ian Cooper
% School of Physics, University of Sydney
% www.physics.usyd.edu.au/teach_res/mp/

clear all
clc
close all


f1 = figure('Color',[0.95 0.9 0.9],'Name','Solving Quadratic Equations', ...
     'NumberTitle','off','Position',[20 50 900 600]);

 t1 = uicontrol(gcf,'Style','text','Position',[250 540 400 30], ...
      'String','Quadratic Equations','FontSize',18, ...
      'HorizontalAlignment','center','FontWeight','bold', ...
      'BackgroundColor',[0.95 0.9 0.9],'ForegroundColor',[0 0.5 1]);
 
 sp1 = subplot('Position',[0.4 0.8 0.1 0.1]); axis off;
 
 t2 = text(0,0,'ax^2 + bx + c = 0','FontName','Arial','FontSize',16, ...
      'FontWeight','bold');
 
% Input Initial Data
a = 1; b = 3; c = -1;
 
Text_a = uicontrol(gcf,'Style','text','Position',[50 380 50 30], ...
          'String','a = ','FontSize',14,'HorizontalAlignment','righ', ...
          'BackgroundColor',[0.95 0.9 0.9]);
      
Edit_a = uicontrol(gcf,'Style','edit','Position',[100 380 100 30], ...
          'String',a,'FontSize',14,'BackgroundColor',[1 1 1], ...
          'Callback','a = str2num(get(Edit_a,''String''));');
      
Text_b = uicontrol(gcf,'Style','text','Position',[50 340 50 30], ...
          'String','b = ','FontSize',14,'HorizontalAlignment','righ', ...
          'BackgroundColor',[0.95 0.9 0.9]);
      
Edit_b = uicontrol(gcf,'Style','edit','Position',[100 340 100 30], ...
          'String',b,'FontSize',14,'BackgroundColor',[1 1 1], ...
          'Callback','b = str2num(get(Edit_b,''String''));');
     
Text_c = uicontrol(gcf,'Style','text','Position',[50 300 50 30], ...
          'String','c = ','FontSize',14,'HorizontalAlignment','righ', ...
          'BackgroundColor',[0.95 0.9 0.9]);
      
Edit_c = uicontrol(gcf,'Style','edit','Position',[100 300 100 30], ...
          'String',c,'FontSize',14,'BackgroundColor',[1 1 1], ...
          'Callback','c = str2num(get(Edit_c,''String''));');      
      
      
 pushbutton_run = uicontrol(gcf,'Style','pushbutton','Position',...
     [100 50 100 40], 'FontSize',12,'FontWeight','bold', 'String','RUN', ...
     'CallBack', 'math_qe_cal');
 
 
 pushbutton_close = uicontrol(gcf,'Style','pushbutton','Position',...
     [300 50 100 40], 'FontSize',12,'FontWeight','bold', 'String','CLOSE', ...
     'CallBack', 'close');
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      