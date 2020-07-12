function beats(action)

% Simulation of beats
% Ian Cooper
% School of Physics, University of Sydney

if (nargin==0)
   action = 'start';
end

switch action


case 'start'
   % ************************************************************************
   % Create plot figure and axes
   % Call local function to create unuser interface controls (uicontrols)
% ************************************************************************    

clc
close all
% clear all

% MAIN FIGURE WINDOW
% INPUTS: location, name ++++++++++++++++++++++++++++++++++++++++++++   
  fig_location =  [2 2 28 18];   
  fig_name     =  'BEATS';   
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig = figure('Name',fig_name,'Visible','off','NumberTitle','off' ...
      ,'Units', 'centimeters','Position', fig_location);

% user inteface controls 
ucH = makeuicontrols;
set(fig,'Backingstore','off','UserData',ucH,'Visible','on',...
        'HandleVisibility','callback');

    
case 'fire'
% *******************************************************************
  % GO!! Get data from Input Boxes
% *******************************************************************      
      
ucH = get(gcf,'UserData');
   
  boxB = get(ucH.boxB,'String');
   if (isempty(boxB))
      boxB = 1100;   % INPUT initial value
   else
      boxB = str2double(boxB);
   end


% *******************************************************************
  % GO!! Do calculations
% ******************************************************************* 

% SETUP 
f(1)  = 1000;              % 1st frequency
f(2)  = boxB;              % 2nd frequency

num = 15001;                % No. of data points for calc. and plotting
                
t_min = 0;                 % time interval values Min x position for image screen
t_max = 2*30/min(f);
t = linspace(t_min,t_max,num);  

%y_1 = zeros(num,1);        % set up matrices for wave function
%y_2 = zeros(num,1);
%y_beat = zeros(num,1);

% CALCULATIONS 
y_1 = sin(2*pi*f(1)*t);
y_2 = sin(2*pi*f(2)*t);
y_beat = y_1 + y_2;  

% *******************************************************************
  % GO!! GRAPHICS: Do the plotting
% ******************************************************************* 

% INPUTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   location_A  =  [8 11 18 6];   % Plot A: location
   location_B  =  [8  2 18 6];   % Plot B: location

   fs = 12;                       % Plots: Font Size
   lineWidth_A1 = 1;              % Plots: Line Widths
   lineWidth_A2 = 1;
   lineWidth_B1 = 2;

   lineColor_A1 = [0 0 0];        % Plots: Line Colors
   lineColor_A2 = [1 0 0];
   lineColor_B1 = [0 0 1];

   title_B = 'BEAT PATTERN';             % Plots: titles - main, axes
   title_Ax = 'time  t (s)';  
   title_Bx = title_Ax;
   title_Ay = 'pressure  P  (a.u)'; 
   title_By = title_Ay;

   x_A1 = t;                       % Plots: data 
   y_A1 = y_1;
   x_A2 = t;
   y_A2 = y_2;

   x_B1 = t;
   y_B1 = y_beat;

  % axis_A = [0 t_max -1.1 1.1];   % Plots: axes limits
  % axis_B = [0 t_max -2.1 2.1];
   axis_A = [0 0.01 -1.1 1.1];   % Plots: axes limits
   axis_B = [0 0.05 -2.1 2.1];
   
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Plot AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
h_plot_A = axes('Units','centimeters','Position',location_A);
set(gcf,'CurrentAxes',h_plot_A);
set(gca,'Fontsize',fs);

plot(x_A1,y_A1,'Color',lineColor_A1,'LineWidth',lineWidth_A1);

hold on
plot(x_A2,y_A2,'Color',lineColor_A2,'LineWidth',lineWidth_A2);

axis(axis_A)
xlabel(title_Ax);
ylabel(title_Ay);

% Plot BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
h_plot_B = axes('Units','centimeters','Position',location_B);
set(gcf,'CurrentAxes',h_plot_B);
set(gca,'Fontsize',fs);

plot(x_B1,y_B1,'Color',lineColor_B1,'LineWidth',lineWidth_B1);

axis(axis_B);
title(title_B);
xlabel(title_Bx);
ylabel(title_By);


 case 'sound_beats'
% *******************************************************************
  % GO!! Get data from Input Boxes
% *******************************************************************      
      
ucH = get(gcf,'UserData');
   
  boxB = get(ucH.boxB,'String');
   if (isempty(boxB))
      boxB = 1100;   % INPUT initial value
   else
      boxB = str2double(boxB);
   end  

% *******************************************************************
  % SOUND Do the calculations
% ******************************************************************* 
   
cf1 = 1000;        
cf2 = boxB;

sf = 22050;                 % sample frequency (Hz)
d = 5.0;                    % duration (s)
n = sf * d;                 % number of samples
s = (1:n) / sf;             % sound data preparation
%s = sin(2 * pi * cf * s);   % sinusoidal modulation
s = sin(2 * pi * cf1 * s)+ sin(2 * pi * cf2 * s);
s = s./max(s);

sound(s, sf);               % sound presentation
pause(d + 0.5);             % waiting for sound end

%filename = 'freq.wav';
%audiowrite(filename,s,sf);

end   %   end case
return;


% *******************************************************************
  % Subsidiary functions
  % Input boxes and Buttons
% ******************************************************************* 

function ucH = makeuicontrols

fcolor = get(gcf,'Color');
fs = 12;     % font size

% BOX AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  
% INPUTS: location, label, unit, initial value ++++++++++++++++++++++
   boxA_label_location  = [1 14 3.5 1];      
   boxA_label = 'f1 = 1000 Hz';   
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

hA_label = uicontrol('Style','text','Units','centimeters',...
           'Position',boxA_label_location,'String',boxA_label,...
           'BackgroundColor',fcolor,'ForegroundColor','k', ...
           'FontSize',fs);

% BOX BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
% INPUTS: location, label, unit, initial value ++++++++++++++++++++++
   boxB_label_location  = [0.2 12 3.5 1];
   boxB_value_location  = [2.6 12.2 2.5 1];
   boxB_unit_location   = [5.0 12 1 1]; 

   boxB_label = 'f2  = ';   
   boxB_value = '1100'; 
   boxB_unit = 'Hz';                    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

h2_label = uicontrol('Style','text','Units','centimeters',...
           'Position',boxB_label_location,'String',boxB_label,...
           'BackgroundColor',fcolor,'ForegroundColor','k', ...
           'FontSize',fs);
      
h2_unit = uicontrol('Style','text','Units','centimeters',...
          'Position',boxB_unit_location,'String',boxB_unit,...
          'BackgroundColor',fcolor,'ForegroundColor','k', ...
          'FontSize',fs);      
 
ucH.boxB = uicontrol('Style','edit','Units','centimeters', ...
           'Position',boxB_value_location, 'String',boxB_value, ...
           'BackgroundColor','white','FontSize',fs);     
         
                 
% BUTTONS -----------------------------------------------------------      
% INPUTS: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   fire_location  = [2 6 3 1.5];     % location of box for GO !! button
   close_location = [2 1 3 1.5];     % location of box for CLOSE
   sound_location = [2 3.5 3 1.5];   % location of box for SOUND

   fire_callback  = 'beats(''fire'')';         % Name of File and case
   sound_callback = 'beats(''sound_beats'')';  % Name of File and caseI
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fire_button = uicontrol('Style','Pushbutton','Units','centimeters',...
              'Position',fire_location,'Callback',fire_callback, ...
              'String','GO ! !','FontSize',fs,'ForegroundColor','b', ...
              'FontWeight','bold');

          
close_button = uicontrol('Style','Pushbutton','Units','centimeters',...
               'Position',close_location,'Callback','close','String', ...
               'Close','FontSize',fs);
           
               
sound_button = uicontrol('Style','Pushbutton','Units','centimeters',...
              'Position',sound_location,'Callback',sound_callback, ...
              'String','SOUND','FontSize',fs);
              
           
return;

