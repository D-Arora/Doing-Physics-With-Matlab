function h = pipeco(action)

if (nargin==0)
   action = 'start';
end

switch action


case 'start',
% ************************************************************************
   % Create plot figure and axes
   % Call local function to create unuser interface controls (uicontrols)
% ************************************************************************    

clc
close all
clear all

% MAIN FIGURE WINDOW
% INPUTS: location, name ++++++++++++++++++++++++++++++++++++++++++++   
  fig_location =  [2 2 28 18];   
  fig_name     =  'Organ Pipe';   
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig = figure('Name',fig_name,'Visible','off','NumberTitle','off' ...
      ,'Units', 'centimeters','Position', fig_location);

% user inteface controls 
ucH = makeuicontrols;

set(fig,'Backingstore','off','UserData',ucH,'Visible','on',...
        'HandleVisibility','callback');

    
case 'fire',
% *******************************************************************
  % GO!! Get data from Input Boxes
% *******************************************************************      

% initial values:
%  mode number, length of pipe, velocity of sound
iv = [2 1 340];

colorbar off
axis off                          % turn off axis for updating plot
%ht = title('xxxxxxxxxxxxxxxxxxxxxxxxxxx');                % remove title when updating plot
%set(ht,'visible', 'off');
hold off

ucH = get(gcf,'UserData');

% mode number
 box1 = get(ucH.box1,'String');
   if (isempty(box1)); box1 = iv(1);   
   else box1 = str2num(box1); end;

% length of pipe
box2 = get(ucH.box2,'String');
if (isempty(box2)); box2 = iv(2);   
else box2 = str2num(box2); end;
 
% velcoity of sound
box3 = get(ucH.box3,'String');
if (isempty(box3)); box3 = iv(3);   
else box3 = str2num(box3); end;   


%   boxB = get(ucH.boxB,'String');
%    if (isempty(boxB))
%       boxB = 1100;   % INPUT initial value
%    else
%       boxB = str2num(boxB);
%    end


% *******************************************************************
  % GO!! Do calculations
% ******************************************************************* 

% SETUP 

global f

N = box1;                 % mode number
L = box2;                 % length of pipe
v = box3;                 % velocity of sound

x = linspace(0,L,100);     % position along pipe
wL = 4*L / (2*N-1);             % wavelength
k = 2*pi/wL;              % propagation constant
f = v / wL;
p = cos(k.*x);            % pressure along pipe


%f(1)  = 1000;              % 1st frequency
%f(2)  = 1100;              % 2nd frequency

num = 5001;                % No. of data points for calc. and plotting
                
t_min = 0;                 % time interval values Min x position for image screen
t_max = 30/min(f);
t = linspace(t_min,t_max,num);  

y_1 = zeros(num,1);        % set up matrices for wave function
y_2 = zeros(num,1);
y_beat = zeros(num,1);

% CALCULATIONS 
%y_1 = sin(2*pi*f(1)*t);
%y_2 = sin(2*pi*f(2)*t);
%y_beat = y_1 + y_2;  

% *******************************************************************
  % GO!! GRAPHICS: Do the plotting
% ******************************************************************* 

% INPUTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   location_A  =  [8 11 18 6];   % Plot A: location
   location_B  =  [8  2 18 6];   % Plot B: location

   fs = 12;                       % Plots: Font Size
             % Plots: Line Widths
   lineWidth_A2 = 1;
   lineWidth_B1 = 2;

   lineColor_A1 = [0 0 0];        % Plots: Line Colors
   lineColor_A2 = [1 0 0];
   lineColor_B1 = [0 0 1];

%    title_B = 'BEAT PATTERN';             % Plots: titles - main, axes
%     
%    %title_Bx = title_Ax;
%    
%    %title_By = title_Ay;
% 
%    
%    x_A2 = x;
%    y_A2 = -p;
% 
%    x_B1 = t;
%    y_B1 = y_beat;
% 
%    axis_A = [0 t_max -1.1 1.1];   % Plots: axes limits
%    axis_B = [0 t_max -2.1 2.1];
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Plot AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
%set(ht,'visible', 'on');
h_plot_A = axes('Units','centimeters','Position',location_A);
set(gcf,'CurrentAxes',h_plot_A);
set(gca,'Fontsize',fs);

x_p = x./L;                       
y_p = p;
     %hp = plot(x_p,y_p);

%lineWidth_p = 2;    
%lineColor_p = [0 0 0];
 %    set(hp,'LineWidth',lineWidth_p);


%hold on

x_p = x./L;                       
y_p = -p;
     %hp = plot(x_p,y_p);
h_area = area(x_p,y_p);
     set(h_area,'FaceColor',[0.2 1 0.41]);

hold on

h_area = area(x_p,-y_p);
     set(h_area,'FaceColor',[0.2 1 0.4]);

%lineWidth_p = 2;    
%lineColor_p = [0 1 0];
 %    set(hp,'Color',lineColor_p,'LineWidth',lineWidth_p);

plot([1 0 0 1],[-1 -1 1 1],'k','LineWidth',10);
    
title_m = 'Organ Pipe: Closed  Open';     
title_x = 'relative postion along pipe'; 
title_y = 'pressure flucatuations'; 

title(title_m);
xlabel(title_x);
ylabel(title_y);

text_1 = '{\itf} = ';
text_2 = num2str(f,'%3.1f');
text_3 = '   Hz';
text_p = [text_1 text_2 text_3];

h_text = text(0.7,1.15,text_p);
     set(h_text,'fontsize',16,'EraseMode','background','BackGroundColor',[1 0 1]);


%global x_p p z xSc, ySc, zSc
% Plot BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
h_plot_B = axes('Units','centimeters','Position',location_B);
set(gcf,'CurrentAxes',h_plot_B);
set(gca,'Fontsize',fs);

%figure(2)


% colormap(gray);
% z = abs(p);
% [xSc, ySc] = meshgrid(x_p',[-1, 1]);
% zSc = [z,z];
% pcolor(xSc,ySc,zSc);
% shading interpret
% axis off

plot([0 0],[0 0]);

colormap(hot);
  [xSc,ySc] = meshgrid(x_p',[-1,+1]);
  zSc = [abs(p);abs(p)];
  pcolor(xSc,ySc,zSc);
  shading interp;
  colorbar
  axis off
%Grayscale(x_p',p')
%axis(axis_B);
%title(title_B);
%xlabel(title_Bx);
%ylabel(title_By);



 case 'sound_beats',
% *******************************************************************
  % GO!! Get data from Input Boxes
% *******************************************************************      
      
% ucH = get(gcf,'UserData');
%    
%   boxB = get(ucH.boxB,'String');
%    if (isempty(boxB))
%       boxB = 1100;   % INPUT initial value
%    else
%       boxB = str2num(boxB);
%    end  

% *******************************************************************
  % SOUND Do the calculations
% ******************************************************************* 
global f

%f = 1000      
sf = 22050;                 % sample frequency (Hz)
d = 5.0;                    % duration (s)
n = sf * d;                 % number of samples
s = (1:n) / sf;             % sound data preparation
s = sin(2 * pi * f * s);   % sinusoidal modulation

s = s./max(s);

sound(s, sf);               % sound presentation
pause(d + 0.5);             % waiting for sound end


end;   %   end case
return;



% *******************************************************************
%   Subsidiary functions
%   Input boxes and Buttons
% ******************************************************************* 

function ucH = makeuicontrols

% Inputs for buttons ++++++++++++++++++++++++++++++++++++++++++++++++
fcolor = get(gcf,'Color');
fs = 12;     % font size


% Inputs for location of input boxes ++++++++++++++++++++++++++++++++
in_num = 3;
box_label_location  = zeros(in_num,4); 
box_value_location  = zeros(in_num,4);
box_unit_location   = zeros(in_num,4);

box_label_location(1,:)  = [1.0  16.0  5.0  1.0]; 
box_value_location(1,:)  = [1.5  15.4  2.5  1.0];
box_unit_location(1,:)   = [4.2  15.2  2.0  1.0]; 

box_label_location(2,:)  = [1.0  14.0  5.0  1.0];
box_value_location(2,:)  = [1.5  13.4  2.5  1.0];
box_unit_location(2,:)  =  [4.2  13.2  2.0  1.0]; 

box_label_location(3,:)  = [1.0  12.0  5.0  1.0];
box_value_location(3,:)  = [1.5  11.4  2.5  1.0];
box_unit_location(3,:)  =  [4.2  11.2  2.0  1.0]; 

% Inputs for label, value, unit +++++++++++++++++++++++++++++++++++++
box_label = {'Mode number, N = 1, 2, 3, ...'; ...
             'Pipe length, L               '; ...
             'velocity of sound, v         '};
box_value = {'1'; '1.2'; '344'};
box_unit =  {'---' ; 'm'; 'm/s'};




% Stepup for input boxes ********************************************
for c = 1:in_num

h_label = uicontrol('Style','text','Units','centimeters',...
           'Position',box_label_location(c,:),'String',char(box_label(c)),...
           'BackgroundColor',fcolor,'ForegroundColor','k', ...
           'FontSize',fs);
       
h_unit = uicontrol('Style','text','Units','centimeters',...
          'Position',box_unit_location(c,:),'String',char(box_unit(c)),...
          'BackgroundColor',fcolor,'ForegroundColor','k', ...
          'FontSize',fs);    
end


ucH.box1 = uicontrol('Style','edit','Units','centimeters', ...
           'Position',box_value_location(1,:), 'String',char(box_value(1)), ...
           'BackgroundColor','white','FontSize',fs); 
       
ucH.box2 = uicontrol('Style','edit','Units','centimeters', ...
           'Position',box_value_location(2,:), 'String',char(box_value(2)), ...
           'BackgroundColor','white','FontSize',fs);  
       
ucH.box3 = uicontrol('Style','edit','Units','centimeters', ...
           'Position',box_value_location(3,:), 'String',char(box_value(3)), ...
           'BackgroundColor','white','FontSize',fs); 



% % BOX AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  
% % INPUTS: location, label, unit, initial value ++++++++++++++++++++++
%    boxA_label_location  = [1 14 3.5 1];      
%    boxA_label = 'f_1 = 1000 Hz';   
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% hA_label = uicontrol('Style','text','Units','centimeters',...
%            'Position',boxA_label_location,'String',boxA_label,...
%            'BackgroundColor',fcolor,'ForegroundColor','k', ...
%            'FontSize',fs);
% 
% % BOX BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
% % INPUTS: location, label, unit, initial value ++++++++++++++++++++++
%    boxB_label_location  = [0.2 12 3.5 1];
%    boxB_value_location  = [2.6 12.2 2.5 1];
%    boxB_unit_location   = [5.0 12 1 1]; 
% 
%    boxB_label = 'f_2  = ';   
%    boxB_value = '1100'; 
%    boxB_unit = 'Hz';                    
% % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
% h2_label = uicontrol('Style','text','Units','centimeters',...
%            'Position',boxB_label_location,'String',boxB_label,...
%            'BackgroundColor',fcolor,'ForegroundColor','k', ...
%            'FontSize',fs);
%       
% h2_unit = uicontrol('Style','text','Units','centimeters',...
%           'Position',boxB_unit_location,'String',boxB_unit,...
%           'BackgroundColor',fcolor,'ForegroundColor','k', ...
%           'FontSize',fs);      
%  
% ucH.boxB = uicontrol('Style','edit','Units','centimeters', ...
%            'Position',boxB_value_location, 'String',boxB_value, ...
%            'BackgroundColor','white','FontSize',fs);     
%          
%                  
% BUTTONS -----------------------------------------------------------      
% INPUTS: +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   fire_location  = [2 6 3 1.5];     % location of box for GO !! button
   close_location = [2 1 3 1.5];     % location of box for CLOSE
   sound_location = [2 3.5 3 1.5];   % location of box for SOUND

   fire_callback  = 'pipe_co(''fire'')';         % Name of File and case
   sound_callback = 'pipe_co(''sound_beats'')';  % Name of File and caseI
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

