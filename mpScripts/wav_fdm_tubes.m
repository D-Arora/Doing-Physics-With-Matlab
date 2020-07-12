function h = air_columns(action)

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
%clear all

global num x dx y p V S v L R R_min R_max flag_inst x_max x_min
global x_trumpet y_trumpet flag_inst Z Y rho_0 R1 R2 L1 L2 x_voice y_voice
global x_oboe y_oboe

num = 2000;                    % number of daya points for plotting
v = 340;                       % speed of sound
rho_0 = 1.2;                   % density of air
x_min = 0; x_max = 0.82;       % length of pipe
x = linspace(x_min,x_max,num);
dx = x(2)-x(1);

p = zeros(1,num);
V = zeros(1,num);
S = zeros(1,num);

% Specify type of instrument ---------------------------------------------
disp('  ');
disp('     Select type of instrument');
disp('  ');
disp('          1:   Cylindical pipe:     Closed   Open');
disp('   ');
disp('          2:   Cylindrical pipe:    Open    Open');
disp('   ');
disp('          3:   Conical pipe (cone): Open     Open');
disp('   ');
disp('          4:   Conical pipe (cone):  Closed   Open');
disp('   ');
disp('          5:   Oboe:                 Closed   Open');
disp('   ');
disp('          6:   Trumpet:              Closed   Open');
disp('   ');
disp('          7:   Organ pipe with hole near end: Open     Open');
disp('   ');
disp('          8:   Human vocal tract "ah ...": two pipe model');
disp('   ');
disp('          9:   Human vocal tract "ah ..."');
disp('   ');


     flag_inst = input('     Enter 1, 2, ... , 9:   ');
disp('   ');

switch flag_inst

    case 1   % Uniform organ pipe: Closed Open
    
   disp('   ');
   L = input('   Length of cylindical pipe (m),  L  =  ');
   disp('   ');
   R = input('   Radius of cylindical pipe (m), R =  ');
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x_max = L + 0.6 * R;
   else
       x_max = L;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1);     
   S = (pi*R^2).* ones(1,num);    
   
   p(1) = 1; V(1) = 0;
   
 case 2   % Uniform organ pipe: Open Open
    
   disp('   ');
   L = input('   Length of cylindical (m),  L  =  ');
   disp('   ');
   R = input('   Radius of cylindical pipe (m), R =  ');
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x_max = L + 0.6 * R;
       x_min = -0.6 * R;
   else
       x_min = 0; x_max = L;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1);     
   S = (pi*R^2).* ones(1,num);    
   
   p(1) = 0; V(1) = 1/S(1);  
   
   
   case 3   % Conical pipe (cone): Open Open
    
   disp('   ');
   L = input('   Length of conical pipe (m),  L  =  ');
   disp('   ');
   R_min = input('   Min radius of conical pipe (m), R_min =  ');
   disp('   ');
   R_max = input('   Max radius of conical pipe (m), R_max =  ');
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x_max = L + 0.6 * R_max; x_min = -0.6 * R_min;
   else
       x_max = L; x_min = 0;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1);
   R = ((R_max-R_min)/L) .*x + R_min;
   S = (pi .* R.^2);   
   
   p(1) = 0; V(1) = 1/S(1);
   
   case 4   % Conical pipe (cone): Closed Open
    
   disp('   ');
   L = input('   Length of conical pipe (m),  L  =  ');
   disp('   ');
   R_min = input('   Min radius of conical pipe (m), R_min =  ');
   disp('   ');
   R_max = input('   Max radius of conical pipe (m), R_max =  ');
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x_max = L + 0.6 * R_max; 
   else
       x_max = L; x_min = 0;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1);
   R = ((R_max-R_min)/L) .*x + R_min;
   S = (pi .* R.^2);   
   
   p(1) = 1; V(1) = 0;
 
   case 5   % Oboe
    load x_oboe;
    load y_oboe;
    x = x_oboe;
    R = y_oboe;
    
    S = pi .* R.^2;
    
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x = x .* (max(x) + 0.6 * max(R))/max(x); 
   end
    
   dx = x(2)-x(1);
   x_min = min(x);
   x_max = max(x);
    
   p(1) = 1; V(1) = 0;  
   
    case 6   % Trumpet
    load x_trumpet;
    load y_trumpet;
    x = x_trumpet;
    R = y_trumpet;
        
    S = pi .* R.^2;
    
   disp('   ');
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x = x .* (max(x) + 0.6 * max(R))/max(x); 
   end
    
   dx = x(2)-x(1);
   x_min = min(x);
   x_max = max(x);
    
   p(1) = 1; V(1) = 0;  
%     
     case 7  % Organ pipe - hole:  open open

     disp('   ');
     disp('   length of organ pipe L  = 0.800 m ');
     disp('   ');
     disp('   radius of organ pipe R  = 0.012 m ');
     disp('   ');
     disp('   hole located 0.160 m from end of organ pipe');
     disp('   ');
     disp('   hole radius  R_hole  = 0.006 m');
     disp('   ');
     
     L = 0.800;  R =  0.012;  
     
   flag_end = input('   Include end correction:  (1) yes   (2) no    ');
   if flag_end ==1
       x_max = L + 0.6 * R;
   else
       x_min = 0; x_max = L;
   end
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1);     
   S = (pi*R^2).* ones(1,num);    
   
   p(1) = 0; V(1) = 1/S(1);   
   
    case 8   % Human vocal tract "ah ...": two pipe model

     disp('   ');
     disp('   length of pipe 1, L1  = 0.090 m ');
     disp('   ');
     disp('   radius of pipe 1, R1  = 5.642e-3 m ');
     disp('   ');
     disp('   length of pipe 2, L2 = 0.080 m');
     disp('   ');
     disp('   radius of pipe 2, R2 = 1.4927e-2 m');
     disp('   ');
     
   L1 = 0.09; R1 = 5.642e-3;
   L2 = 0.08; R2 = 1.4927e-2;  
   L = L1+L2; R = R1;
   
   
   x_min = 0; x_max = L;
   x = linspace(x_min,x_max,num);
   dx = x(2)-x(1);     
   S = (pi*R^2).* ones(1,num);    
   for c = 1 : num                 % can rewrite using logical functions
   if x(c) > L1
       R = R2;
       S(c) = (pi*R^2);   
   end
   end
   
   p(1) = 1; V(1) = 0;   
   
   case 9   % Human voice tract: "ah ..."
    load x_voice;
    load y_voice;
    x = x_voice;
    R = y_voice;
    S = pi .* R.^2;
    
   disp('   ');
   x_min = 0; x_max = max(x);
    
   dx = x(2)-x(1);
   x_min = min(x);
   x_max = max(x);
    
   p(1) = 1; V(1) = 0;  
   
end
    

% MAIN FIGURE WINDOW
% INPUTS: location, name ++++++++++++++++++++++++++++++++++++++++++++   
  fig_location =  [2 2 28 18];   
  fig_name     =  'Air Columns';   
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

axis off

ucH = get(gcf,'UserData');
   
  boxB = get(ucH.boxB,'String');
   if (isempty(boxB))
      boxB = 100;   % INPUT initial value
   else
      boxB = str2num(boxB);
   end


% *******************************************************************
  % GO!! Do calculations
% ******************************************************************* 

% SETUP 
f  = boxB;              % 2nd frequency

global num x y dx p V S v L R R_min R_max flag_inst x_max x_min
global x_trumpet y_trumpet flag_inst Z Y rho_0 R1 R2 L1 L2 x_voice y_voice
global x_oboe y_oboe

k = 2*pi*f/v;

% Impedance Z and Admittance Y --------------------------------------------------

Z = (2*pi*f*rho_0*dx) ./ S;
Y = (2*pi*f*dx/(rho_0*v^2)) .* S;
R_hole = 0.006;
Y_hole = R_hole /(2.80*rho_0*f);

for c = 2:num
   phalf = p(c-1) + 0.5 * V(c-1) * Z(c-1) ;
   Vhalf = V(c-1) - 0.5 * p(c-1) * Y(c-1);
   
   if flag_inst==6
       c1 = 1599; c2 = 1616;
       if (c > c1) & (c < c2)
          Vhalf = Vhalf + p(c-1)*Y_hole;
       end
   end

   p(c) = p(c-1) + Vhalf * Z(c-1);
   V(c) = V(c-1) - phalf * Y(c-1);
   
   if flag_inst==6
       c1 = 1599; c2 = 1616;
       if (c > c1) & (c < c2)
          V(c) = V(c) + p(c)*Y_hole;
       end
   end   
end
p_max = max(abs(p));
p = p./p_max;

% *******************************************************************
  % GO!! GRAPHICS: Do the plotting
% ******************************************************************* 

% INPUTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   location_A  =  [8 11 18 6];   % Plot A: location
   location_B  =  [8  2 18 6];   % Plot B: location

   fs = 12;                       % Plots: Font Size
 

% Plot AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
% Shape of pipe
   h_plot_A = axes('Units','centimeters','Position',location_A);
   set(gcf,'CurrentAxes',h_plot_A);
   set(gca,'Fontsize',fs);

lineWidth_p = 8;
lineColor_p = [0 0 0];

title_p = 'Shape of pipe';             

switch flag_inst
    case 1
    x_p = [ L  0 0 L];
    y_p = [-R -R R R];
    plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    set(gca,'xLim',[0 x_max*1.2]);
    
    case 2
    x_p = [  0   L    0  L];
    y_p = [ -R  -R    R  R];
    plot(x_p(1:2),y_p(1:2),'Color',lineColor_p,'LineWidth',lineWidth_p);
    hold on
    plot(x_p(3:4),y_p(3:4),'Color',lineColor_p,'LineWidth',lineWidth_p);
    set(gca,'xLim',[1.2*x_min x_max*1.2]);
    
    case 3
    x_p = [0       L       0      L];
    y_p = [-R_min  -R_max  R_min  R_max];
    plot(x_p(1:2),y_p(1:2),'Color',lineColor_p,'LineWidth',lineWidth_p);
    hold on
    plot(x_p(3:4),y_p(3:4),'Color',lineColor_p,'LineWidth',lineWidth_p);
    set(gca,'xLim',[1.2*x_min x_max*1.2]);
    
    case 4
    x_p = [L       0       0      L];
    y_p = [-R_max  -R_min  R_min  R_max];
    plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    %hold on
    %plot(x_p(3:4),y_p(3:4),'Color',lineColor_p,'LineWidth',lineWidth_p);
    set(gca,'xLim',[1.2*x_min x_max*1.2]);
    
    case 5
    lineWidth_p = 3;    
    x_p = x_oboe;
    y_p = y_oboe;
    plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    hold on
    plot(x_p,-y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    
    x_p = [0 0];
    y_p = [-y_oboe(1) y_oboe(1)];
    plot(x_p,-y_p,'Color',lineColor_p,'LineWidth',2*lineWidth_p);
    set(gca,'xLim',[1.2*x_min x_max*1.2]); 
    
    case 6
    lineWidth_p = 3;    
    x_p = x_trumpet;
    y_p = y_trumpet;
    plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    hold on
    plot(x_p,-y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    
    x_p = [0 0];
    y_p = [-y_trumpet(1) y_trumpet(1)];
    plot(x_p,-y_p,'Color',lineColor_p,'LineWidth',2*lineWidth_p);
    set(gca,'xLim',[1.2*x_min x_max*1.2]); 
    
    case 7
    x_p = [ 0 L]; y_p = [-R -R];
       plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    hold on
    x_p = [ 0 x(c1)]; y_p = [R R];
       plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    x_p = [ x(c2) L]; y_p = [R R];
       plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    
    set(gca,'xLim',[0 x_max*1.2]);
    
    case 8
    x_p = [ L1+L2  L1  L1  0  0  L1 L1 L1+L2];
    y_p = [-R2    -R2 -R1 -R1 R1 R1 R2 R2 ];
    plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    set(gca,'xLim',[0 x_max*1.2]);    
    
    case 9
    lineWidth_p = 3;    
    x_p = x_voice;
    y_p = y_voice;
    plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    hold on
    plot(x_p,-y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
    
    x_p = [0 0];
    y_p = [-y_voice(1) y_voice(1)];
    plot(x_p,-y_p,'Color',lineColor_p,'LineWidth',2*lineWidth_p);
    set(gca,'xLim',[1.2*x_min x_max*1.2]); 
end

title(title_p);
   

% Plot BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
% Pressure distribution along tube

h_plot_B = axes('Units','centimeters','Position',location_B);
set(gcf,'CurrentAxes',h_plot_B);
set(gca,'Fontsize',fs);

lineWidth_p = 2;
lineColor_p = [0 0 1];

title_1 = ('end pressure at (+)  =  ');
title_2 = num2str(p(end),'%0.3f');
title_p = [title_1 title_2];            
title_y = 'pressure  P  (a.u)'; 
title_x = 'position  x  (m)';

plot(x_max,0,'+','MarkerEdgeColor',[1 0 0],'Markersize',12)
hold on

x_p = x;
y_p = p;

plot(x_p,y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);

lineColor_p = [1 0 0];
plot(x_p,-y_p,'Color',lineColor_p,'LineWidth',lineWidth_p);
set(gca,'xLim',[x_min*1.2 x_max*1.2]);
%set(gca,'yLim',[-1.1 1.1]);

%axis(axis_B);
h_title = title(title_p);
set(h_title,'BackgroundColor',[0.8 0.8 0.8],'EraseMode','background');
xlabel(title_x);
ylabel(title_y);


 case 'sound_beats',
% *******************************************************************
  % GO!! Get data from Input Boxes
% *******************************************************************      
      
ucH = get(gcf,'UserData');
   
  boxB = get(ucH.boxB,'String');
   if (isempty(boxB))
      boxB = 1100;   % INPUT initial value
   else
      boxB = str2num(boxB);
   end  

% *******************************************************************
  % SOUND Do the calculations
% ******************************************************************* 
   
%cf1 = 1000;        
cf2 = boxB;
%cf1
%cf2
sf = 22050;                 % sample frequency (Hz)
d = 5.0;                    % duration (s)
n = sf * d;                 % number of samples
s = (1:n) / sf;             % sound data preparation
s = sin(2 * pi * cf * s);   % sinusoidal modulation
%s = sin(2 * pi * cf1 * s)+ sin(2 * pi * cf2 * s);
s = s./max(s);

sound(s, sf);               % sound presentation
pause(d + 0.5);             % waiting for sound end
end;   %   end case
return;


% *******************************************************************
  % Subsidiary functions
  % Input boxes and Buttons
% ******************************************************************* 

function ucH = makeuicontrols

fcolor = get(gcf,'Color');
fs = 12;     % font size

% BOX BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
% INPUTS: location, label, unit, initial value ++++++++++++++++++++++
   boxB_label_location  = [0.2 9.0 3.5 1];
   boxB_value_location  = [2.6 9.2 2.5 1];
   boxB_unit_location   = [5.0 9.0 1 1]; 

   boxB_label = 'f  = ';   
   boxB_value = '100'; 
   boxB_unit = 'Hz';                    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

h2_label = uicontrol('Style','text','Units','centimeters',...
           'Position',boxB_label_location,'String',boxB_label,...
           'BackgroundColor',fcolor,'ForegroundColor','k', ...
           'FontSize',16);
      
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

   fire_callback  = 'air_columns(''fire'')';         % Name of File and case
   sound_callback = 'air_columns(''sound_beats'')';  % Name of File and caseI
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

