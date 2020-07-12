function satellite(action)
clc
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
clear 

% MAIN FIGURE WINDOW
% INPUTS: location, name ++++++++++++++++++++++++++++++++++++++++++++   
  fig_location =  [2.0 2.0 28.0 18.5];   
  fig_name     =  'ORBITS ';   
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig = figure('Name',fig_name,'Visible','off','NumberTitle','off' ...
      ,'Units', 'centimeters','Position', fig_location);

% user inteface controls 
ucH = makeuicontrols;
     set(fig,'Backingstore','off','UserData',ucH,'Visible','on',...
        'HandleVisibility','callback');

    
case 'fire'
% *******************************************************************
  % Get data from Input Boxes
% *******************************************************************  

iv= [8000 1.2 0 90 5000 1000];
axis off

ucH = get(gcf,'UserData');

 box1 = get(ucH.box1,'String');
   if (isempty(box1))
       box1 = iv(1);   
   else
       box1 = str2double(box1);
   end
  
box2 = get(ucH.box2,'String');
  if (isempty(box2))
    box2 = iv(2);   
  else
      box2 = str2double(box2);
  end
   
box3 = get(ucH.box3,'String');
  if (isempty(box3))
    box3 = iv(3);   
  else
      box3 = str2double(box3);
  end   

box4 = get(ucH.box4,'String');
  if (isempty(box4))
    box4 = iv(4);   
  else
      box4 = str2double(box4);
  end  

box5 = get(ucH.box5,'String');
  if (isempty(box5))
    box5 = iv(5);   
  else
      box5 = str2double(box5);
  end  

box6 = get(ucH.box6,'String');
  if (isempty(box6))
    box6 = iv(6);   
  else
      box6 = str2double(box6);
  end  
   

% *******************************************************************
  % START CALCULATIONS
% ******************************************************************* 

% Input and Constants -----------------------------------------------
v_0 = box1;                % initial launch velocity (m/s)
R_0 = box2;                % Initial displacement of satellite / R_E
lat = box3;                % latitude for launch of satellite (deg)
angle_launch = box4;       % launch angle of satellite (deg) 
t_max = box5;              % max simulation time (s)
nt = box6;                 % number of time steps

R_E = 6.38e6;              % Earth's radius
M_E = 5.98e24;             % Earth's mass
G = 6.67e-11;              % Universal gravitation constant

%  Setup and initialize variables -----------------------------------
v_0x = v_0 * cosd(angle_launch);         % initial velcoity: x cpt.
v_0y = v_0 * sind(angle_launch);         % initial velcoity: y cpt.
R_0x = R_E * R_0 * cosd(lat);            % initial displacement: x cpt.
R_0y = R_E * R_0 * sind(lat);            % initial displacement: y cpt.

t = linspace(0,t_max,nt);                % time
dt = t(2) - t(1);                        % time step 

v = zeros(1,nt);


% ********************************************************************
%  START NUMERICAL CALCULATIONS
% ********************************************************************

% time steps 1 and 2 (time 1 = 0) -----------------------------------
R(1) = R_0x + 1i * R_0y;
v(1) = v_0x + 1i * v_0y;
a(1) = -(G*M_E / abs(R(1))^3) * R(1);
v(2) = v(1) + a(1) * dt;
R(2) = R(1) + v(1) * dt + 0.5 * a(1) * dt^2;


% time steps 3, 4, 5, .... to calculate position --------------------
for c = 3 : nt
   R(c) = -G*M_E*dt^2*R(c-1)/abs(R(c-1))^3 + 2*R(c-1) - R(c-2);
end

Rx = real(R); Ry = imag(R);

x_max = max(real(R))/R_E; x_min = min(real(R))/R_E;   % max x    min x
y_max = max(imag(R))/R_E; y_min = min(imag(R))/R_E;   % max y    min y
a = (x_max-x_min)/2;                                  % semimajor axis
b = (y_max-y_min)/2;                                  % semiminor axis            
e = abs((x_max+x_min)/(2*a));                         % eccentricity

% velocity ----------------------------------------------------------
vx = zeros(1,nt); vy = zeros(1,nt);
vx(1) = (Rx(2)-Rx(1))/dt; vy(1) = (Ry(2)-Ry(1))/dt; 
vx(nt) = (Rx(end)-Rx(end-1))/dt; vy(nt) = (Ry(end)-Ry(end-1))/dt; 
for c = 2: nt-1
    vx(c) = (Rx(c+1)-Rx(c-1))/(2*dt);
    vy(c) = (Ry(c+1)-Ry(c-1))/(2*dt);
end
    v = sqrt(vx.^2+ vy.^2);
    
% acceleration ------------------------------------------------------
ax = zeros(1,nt); ay = zeros(1,nt);
ax(1) = (vx(2)-vx(1))/dt; ay(1) = (vy(2)-vy(1))/dt; 
ax(nt) = (vx(end)-vx(end-1))/dt; ay(nt) = (vy(end)-vy(end-1))/dt; 
for c = 2: nt-1
    ax(c) = (vx(c+1)-vx(c-1))/(2*dt);
    ay(c) = (vy(c+1)-vy(c-1))/(2*dt);
end
    acc = sqrt(ax.^2+ ay.^2);   


% energy / mass -----------------------------------------------------
K_m = 0.5*v.^2; U_Gm = -(G*M_E)./abs(R); E_m = K_m + U_Gm;  

L_m = v.*abs(R);              % angular momentum / mass

% circular orbit ----------------------------------------------------
v_orb = sqrt(G*M_E/(a*R_E));        % orbital velocity
g = G*M_E/(a*R_E)^2;                % acceleratin due to gravity
T = 2*pi*(a*R_E)/v_orb;             % period
K_mA = 0.5 * v(end)^2;              % kinetic energy / mass
U_GmA = -(G*M_E)/(a*R_E);           % potenetial energy / mass 
E_mA = K_mA + U_GmA;                % total energy 


% Output parameters at end of simulation --------------------------------
x_end    = real(R(end))/R_E;          % x coordinate of position
y_end    = imag(R(end))/R_E;          % y coordinate of position
R_end    = abs(R(end))/R_E;           % distance from Earth
R_angle  = rad2deg(angle(R(end)));    % direction of displacement vector

v_end = v(end);                        % velocity magnitude
v_angle = atan2d(vy(end),vx(end));     % velocity direction
vx_end   = vx(end);                    % velocity components
vy_end   = vy(end);

a_end = acc(end);                      % acceleration magnitude
a_angle = atan2d(ay(end),ax(end));     % acceleration direction
ax_end   = ax(end);                    % acceleration components
ay_end   = ay(end);


Km_end  = K_m(end);                   % kinetic energy
UGm_end = U_Gm(end);                  % gravitational potential energy
Em_end  = E_m(end);                   % total energy


% *******************************************************************
%  GRAPHICS
% ******************************************************************* 

 % Plot SETUP -------------------------------------------------------
figure(1)
% pos = [0.05 0.05 0.55 0.7];
%   set(gcf,'Units','normalized');
%   set(gcf,'Position',pos);
%   set(gcf,'color','w'); 
 location_A  =  [8 3 18 14];   % Plot A: location
     h_plot_A = axes('Units','centimeters','Position',location_A);

set(gcf,'CurrentAxes',h_plot_A);

fs = 12;                      % font size
set(gca,'Fontsize',fs);

limit_p = ceil(max(abs(R))/R_E);
step_p = 0.5;
if limit_p > 4.9 ; step_p = 1; end
if limit_p > 10.0; step_p = 2; end
if limit_p > 20.0; step_p = 5; end
if limit_p > 50.0; step_p = 10; end
     set(gca,'Xlim',[-limit_p limit_p]);
     set(gca,'Xtick',-limit_p:step_p:limit_p);
     set(gca,'Ylim',[-limit_p limit_p]);
     set(gca,'Ytick',-limit_p:step_p:limit_p);

% Plot: x, y trajectory ---------------------------------------------
%set(gcf,'Name','Orbit');


x_p = real(R)./R_E;     % x data for plot
y_p = imag(R)./R_E;     % y data for plot
     H = plot(x_p,y_p);

lineWidth_p = 1.5;
     set(H,'lineWidth',lineWidth_p);

title_x = 'x position  (m)';
title_y = 'y position  (m)';
     xlabel(title_x);
     ylabel(title_y);
     
hold on

color_p = [0 0 1];
color_p1 = [1 1 1];
size_p = 4;
ns = nt /50;
h_plot_A2 = plot(x_p(1:ns:end),y_p(1:ns:end),'o');
     set(h_plot_A2,'MarkerEdgeColor',color_p,'MarkerFaceColor',color_p1, ...
        'MarkerSize',size_p);

color_p = [1 0 0];
size_p = 8;    
h_plot_A2 = plot(x_p(end),y_p(end),'o');
     set(h_plot_A2,'MarkerEdgeColor',color_p,'MarkerFaceColor',color_p, ...
        'MarkerSize',size_p);    

% plot Earth --------------------------------------------------------
%axis_dim = [-4 4 -4 4];
color_face = [1 0.6 0.3];
h_rectangle = rectangle('Position', [-1,-1,2,2],'Curvature',[1 1]);
set(h_rectangle,'FaceColor',color_face,'EdgeColor',color_face);

% plot points for focus and center of orbit -------------------------
color_p = [0 0 1];
color_p1 = [0 0 1];
size_p = 10;
plot([0 -2*e*a],[0 0],'+','MarkerEdgeColor',color_p,...
          'MarkerFaceColor',color_p1,'Markersize',size_p);
 
color_p = [1 0.6 0.3]; color_p1 = [0 0 0];          
plot(-e*a,0,'o','MarkerEdgeColor',color_p,...
          'MarkerFaceColor',color_p1,'Markersize',size_p);          

axis equal
box on
grid on
set(gca,'Fontsize',fs);

% Plot: velocity and orbit radius / time ----------------------------
figure(3);
  pos = [0.6 0.65 0.25 0.22];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
%set(gcf,'Name','Velocity','Number','off');

x_p = t;                % x data for plot
y_p1 = v;               % y1 data for plot
y_p2 = abs(R)./R_E;     % y2 data for plot
     [AX, H1, H2] = plotyy(x_p,y_p1,x_p,y_p2);

title_x = 'time  t  (s)';
title_y1 = 'velocity  v  (m/s)';
title_y2 = 'orbit radius  R  (m)';
lineWidth_p = 2;
color_p1 = [0 0 0]; color_p2 = [0 0 1];
     xlabel(title_x);
     set(get(AX(1),'Ylabel'),'String',title_y1);
     set(get(AX(2),'Ylabel'),'String',title_y2);
     set(AX(1),'YColor', color_p1);
     set(AX(2),'YColor', color_p2);
     set(H1,'LineWidth',lineWidth_p, 'Color',color_p1);
     set(H2,'LineWidth',lineWidth_p, 'Color',color_p2);


% setting the ytick marks for both Y axes     
limit_ymax = 10^(ceil(log10(max(v))));
if limit_ymax > 30000; limit_ymax = 30000; end
limit_ymin = 10^(floor(log10(min(v))));
limit_dy = (limit_ymax - limit_ymin)/9;
ytick_p = limit_ymin:limit_dy:limit_ymax;
     set(AX(1),'Ytick',ytick_p);
     
limit_ymax = 10^(ceil(log10(max(abs(R)./R_E))));
limit_ymin = 10^(floor(log10(min(abs(R)./R_E))));
limit_dy = (limit_ymax - limit_ymin)/9;
ytick_p = limit_ymin:limit_dy:limit_ymax;
    set(AX(2),'YLim',[limit_ymin limit_ymax]);
     set(AX(2),'Ytick',ytick_p);     
     
% Plot: angular momentum / time -------------------------------------

figure(4)
  pos = [0.6 0.35 0.25 0.22];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
%set(gcf,'Name','Angular Momentum','Number','off');

x_p = t;     % x data for plot
y_p = L_m;     % y data for plot
     H = plot(x_p,y_p);

lineWidth_p = 2;     
     set(H,'lineWidth',lineWidth_p);
     
title_x = 'time  t  (s)';
title_y = 'Angular Momentum L/m  (m  ^2/s)';
     xlabel(title_x);
     ylabel(title_y);

limit_y = [0 2*max(L_m)] ;    
     set(gca,'Ylim',limit_y);

% Plot: energy / time -----------------------------------------------
figure(5)
  pos = [0.6 0.05 0.25 0.22];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w'); 
%set(gcf,'Name','Energy','Number','off');
hold off

x_p = t;     % x data for plot
y_p = K_m;     % y data for plot
     H = plot(x_p,y_p);

lineWidth_p = 2;
     set(H,'lineWidth',lineWidth_p);

title_x = 'time  t  (s)';
title_y = 'Energy / mass  (J.kg  {^-1})';
     xlabel(title_x);
     ylabel(title_y);

hold on

x_p = t;     % x data for plot
y_p = U_Gm;     % y data for plot
     H = plot(x_p,y_p);
     
color_p = [0 0 0];
     set(H,'Color',color_p,'lineWidth',lineWidth_p)

x_p = t;       % x data for plot
y_p = E_m;     % y data for plot
     H = plot(x_p,y_p);
     
color_p = [1 0 1];
     set(H,'Color',color_p,'lineWidth',lineWidth_p);
     
legend('K','U_{Gm}','E_m');


% *******************************************************************
%   OUTPUT RESULTS TO SCREEN
% ******************************************************************* 
%clc
disp('  ');
disp('Input Parameters');
fprintf('   launch velocity, v_0    =  %0.1f  m/s  \n',v_0);
fprintf('   launch radius, R_0      =  %0.2f  ---  \n',R_0);
fprintf('   launch latitude         =  %0.1f  deg  \n',lat);
fprintf('   launch angle            =  %0.1f deg  \n',angle_launch);
fprintf('   simulation time, tamx   =  %0.1f  s  \n',t_max);
disp('   ')
disp('Output Parameters at end of simulation');
fprintf('Displacement vector: magnitude R =  %0.2f   \n',R_end);
fprintf('Displacement vector: direction theta_R =  %0.2f  deg  \n',R_angle);
fprintf('   R_x =  %0.2f   \n',x_end); 
fprintf('   R_y =  %0.2f   \n',y_end);
disp('  ')
fprintf('Velocity vector: magnitude v =  %0.2f  m/s  \n',v_end);
fprintf('Velocity vector: direction theta_v =  %0.2f  deg  \n',v_angle);
fprintf('   v_x =  %0.2f  m/s  \n',vx_end); 
fprintf('   v_y =  %0.2f  m/s  \n',vy_end);

disp('  ')
fprintf('Acceleration vector: magnitude a =  %0.2f  m/s^2  \n',a_end);
fprintf('Acceleration vector: direction theta_a =  %0.2f  deg  \n',a_angle);
fprintf('   a_x =  %0.2f  m/s^2  \n',ax_end); 
fprintf('   a_y =  %0.2f  m/s^2  \n',ay_end);
disp('  ')
fprintf('KE/m =  %2.2e  J/kg  \n',Km_end); 
fprintf('UG/m =  %2.2e  J/kg  \n',UGm_end);
fprintf('E/m  =  %2.2e  J/kg  \n',Em_end);
end   %   end case

return;


% *******************************************************************
%   Subsidiary functions
%   Input boxes and Buttons
% ******************************************************************* 

function ucH = makeuicontrols

% Inputs for buttons ++++++++++++++++++++++++++++++++++++++++++++++++
fcolor = get(gcf,'Color');
fs = 12;     % font size
fire_location  = [1.5 3.0 3.0 1.5];     % location of box for GO !! button
close_location = [1.5 1.0 3.0 1.5];     % location of box for CLOSE
fire_callback  = 'mec_satellite_gui(''fire'')';         % Name of File

% Inputs for location of input boxes ++++++++++++++++++++++++++++++++
in_num = 6;
box_label_location  = zeros(in_num,4); 
box_value_location  = zeros(in_num,4);
box_unit_location   = zeros(in_num,4);

box_label_location(1,:)  = [1.0  17.0  4.0  1.0]; 
box_value_location(1,:)  = [1.5  16.4  2.5  1.0];
box_unit_location(1,:)   = [4.2  16.2  2.0  1.0]; 

box_label_location(2,:)  = [1.0  15.0  5.0  1.0];
box_value_location(2,:)  = [1.5  14.4  2.5  1.0];
box_unit_location(2,:)  =  [4.2  14.2  2.0  1.0]; 

box_label_location(3,:)  = [1.0  13.0  2.6  1.0];
box_value_location(3,:)  = [1.5  12.4  2.5  1.0];
box_unit_location(3,:)  =  [4.2  12.2  2.0  1.0]; 

box_label_location(4,:)  = [1.0  11.0  3.5  1.0];
box_value_location(4,:)  = [1.5  10.4  2.5  1.0];
box_unit_location(4,:)  =  [4.2  10.2  2.0  1.0]; 

box_label_location(5,:)  = [1.0  9.0  2.9  1.0];
box_value_location(5,:)  = [1.5  8.4  2.5  1.0];
box_unit_location(5,:)  =  [4.2  8.2  2.0  1.0]; 

box_label_location(6,:)  = [1.0  7.0  3.3  1.0];
box_value_location(6,:)  = [1.5  6.4  2.5  1.0];
box_unit_location(6,:)  =  [4.2  6.2  2.0  1.0]; 
  
% Inputs for label, value, unit +++++++++++++++++++++++++++++++++++++
box_label = {'Launch Velocity'; 'Launch Radius / R_E'; 'Latitude';...
             'Launch Angle';'max time'; '# time steps'};
box_value = {'8000'; '1.2'; '0'; '90'; '5000'; '1000'};
box_unit =  {'m/s' ; '---'; 'deg';'deg'; 's'; '---'};

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

ucH.box4 = uicontrol('Style','edit','Units','centimeters', ...
           'Position',box_value_location(4,:), 'String',char(box_value(4)), ...
           'BackgroundColor','white','FontSize',fs);
       
ucH.box5 = uicontrol('Style','edit','Units','centimeters', ...
           'Position',box_value_location(5,:), 'String',char(box_value(5)), ...
           'BackgroundColor','white','FontSize',fs);
       
ucH.box6 = uicontrol('Style','edit','Units','centimeters', ...
           'Position',box_value_location(6,:), 'String',char(box_value(6)), ...
           'BackgroundColor','white','FontSize',fs);
                 
% BUTTONS -----------------------------------------------------------      
fire_button = uicontrol('Style','Pushbutton','Units','centimeters',...
              'Position',fire_location,'Callback',fire_callback, ...
              'String','GO ! !','FontSize',fs,'ForegroundColor','b', ...
              'FontWeight','bold');

          
close_button = uicontrol('Style','Pushbutton','Units','centimeters',...
               'Position',close_location,'Callback','close all','String', ...
               'Close','FontSize',fs);
           
return;

