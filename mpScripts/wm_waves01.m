function h = wm_waves01(action)

if (nargin==0)
   action = 'start';
end

switch action

case 'start'

% Create plot figure and axes.
fig = figure('Name','Travelling Sine Wave','Visible','off','NumberTitle','off');
axes('Units','normalized','Position',[0.10 0.30 0.85 0.65]);

% Call local function to create uicontrols.
ucH = makeuicontrols;

set(fig,'Backingstore','off','UserData',ucH,'Visible','on',...
        'HandleVisibility','callback');

% Finished with polytool initialization.

case 'edittext'

case 'fire'
   ucH = get(gcf,'UserData');
 
   box1 = get(ucH.box1,'String');
   if (isempty(box1))
      box1 = 10;
   else
      box1 = str2double(box1);
   end

   box2 = get(ucH.box2,'String');
   if (isempty(box2))
      box2 = 100;
   else
      box2 = str2double(box2);
  end
  
  box3 = get(ucH.box3,'String');
   if (isempty(box3))
      box3 = 100;
   else
      box3 = str2double(box3);
  end
  
  box4 = get(ucH.box4,'String');
   if (isempty(box4))
      box4 = 0;
   else
      box4 = str2double(box4);
  end
  
  box5 = get(ucH.box5,'String');
   if (isempty(box5))
      box5 = 0;
   else
      box5 = str2double(box5);
   end

   box6 = get(ucH.box6,'String');
   if (isempty(box6))
      box6 = 1;
   else
      box6 = str2double(box6);
   end
   
   
   
% calculate sine function -----------------------------------------------   
n = 400;
xmin = 0;
xmax = 100;
ymin = -10;
ymax = 10;

A = box1;    % amplitude
wl = box2;   % wavelength
T = box3;   % period
phi = pi*box4;   % initial phase angle
time = box5;  % time evolved

x = linspace(xmin,xmax,n);

if box6 == 1
y = A.*sin(2*pi.*(x./wl - time/T) + phi);
x1 = x(n/2)+(wl/T)*time;
y1 = A.*sin(2*pi.*(x1/wl - time/T) + phi);
else
y = A.*sin(2*pi.*(x./wl + time/T) + phi);    
x1 = x(n/2)-(wl/T)*time;
y1 = A.*sin(2*pi.*(x1/wl + time/T) + phi); 
end

% Graphics ------------------------------------------------------------   
plot(x,y,'b','LineWidth',2)
axis([0 xmax ymin ymax])
%    title(['Range: ',num2str(range),' ft   Flight time: ',num2str(ttime),' s'])
title('y = A sin(2\pi(x/\lambda -/+ t/T) + \phi)') 
xlabel('position  x')
ylabel('wavefunction  y')
grid on
hold on
plot(x1,y1,'ro','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8)
hold off
FS = 14;
set(gca,'fontsize',FS)
end

function ucH = makeuicontrols
% MAKEUICONTROLS   Local function to create uicontrols for polytool. 
    pos = [0.05 0.05 0.6 0.7];
    set(gcf,'Units','normalized');
    set(gcf,'Position',pos);
    FS = 14;
box on
%offset = 0;
fcolor = get(gcf,'Color');
dx = 0.10;
 
% box1------------------------------------------------------------------------
field = [0.08 0.15 0.10 0.05];

uicontrol('Style','text','Units','normalized',...
      'Position',field,'String','A: 0 to 10',...
	  'BackgroundColor',fcolor,'ForegroundColor','b', ...
      'fontsize',FS);

field(1) = field(1) + dx;

ucH.box1 = uicontrol('Style','edit','Units','normalized','Position',field,...
         'String','10','BackgroundColor','white','UserData',100,...
         'CallBack','wm_waves01(''edittext'')');

set(ucH.box1,'Fontsize',FS);   

% box2------------------------------------------------------------------------
field = [0.08 0.05 0.10 0.06];

uicontrol('Style','text','Units','normalized',...
      'Position',field,'String','wL: 1 to 100',...
	  'BackgroundColor',fcolor,'ForegroundColor','b', ...
      'fontsize',FS);

field(1) = field(1) + dx;

ucH.box2 = uicontrol('Style','edit','Units','normalized','Position',field,...
         'String','100','BackgroundColor','white','UserData',100,...
         'CallBack','wm_waves01(''edittext'')');
set(ucH.box2,'Fontsize',FS);   

% box3------------------------------------------------------------------------
field = [0.30 0.15 0.10 0.05];

uicontrol('Style','text','Units','normalized',...
      'Position',field,'String','T: 1 to 100',...
	  'BackgroundColor',fcolor,'ForegroundColor','b', ...
      'fontsize',FS);

field(1) = field(1) + dx;

ucH.box3 = uicontrol('Style','edit','Units','normalized','Position',field,...
         'String','100','BackgroundColor','white','UserData',100,...
         'CallBack','wm_waves01(''edittext'')');
set(ucH.box3,'Fontsize',FS);        
     
% box4------------------------------------------------------------------------
field = [0.30 0.05 0.10 0.06];

uicontrol('Style','text','Units','normalized',...
      'Position',field,'String','phi/pi: 0 to 2',...
	  'BackgroundColor',fcolor,'ForegroundColor','b', ...
      'fontsize',FS);

field(1) = field(1) + dx;

ucH.box4 = uicontrol('Style','edit','Units','normalized','Position',field,...
         'String','0','BackgroundColor','white','UserData',0,...
         'CallBack','wm_waves01(''edittext'')');
set(ucH.box4,'Fontsize',FS);   
     
% box5------------------------------------------------------------------------
field = [0.5 0.15 0.10 0.05];

uicontrol('Style','text','Units','normalized',...
      'Position',field,'String','time t',...
	  'BackgroundColor',fcolor,'ForegroundColor','b', ...
      'fontsize',FS);

  
field(1) = field(1) + dx;

ucH.box5 = uicontrol('Style','edit','Units','normalized','Position',field,...
         'String','0','BackgroundColor','white','UserData',0,...
         'CallBack','wm_waves01(''edittext'')');
set(ucH.box5,'Fontsize',FS);   

% box6------------------------------------------------------------------------
field = [0.5 0.05 0.10 0.06];
% 
uicontrol('Style','text','Units','normalized',...
       'Position',field,'String','R: 1 or L: 2',...
 	  'BackgroundColor',fcolor,'ForegroundColor','b', ...
      'fontsize',FS);
 
field(1) = field(1) + dx;
 
ucH.box6 = uicontrol('Style','edit','Units','normalized','Position',field,...
          'String','1','BackgroundColor','white','UserData',0,...
          'CallBack','wm_waves01(''edittext'')');
set(ucH.box6,'Fontsize',FS);   

% buttons------------------------------------------------------------------------     
     
fire_button = uicontrol('Style','Pushbutton','Units','normalized',...
               'Position',[0.75 0.15 0.12 0.05],'Callback','wm_waves01(''fire'')','String','GO ! !');
set(fire_button,'Fontsize',FS);    
           
close_button = uicontrol('Style','Pushbutton','Units','normalized',...
               'Position',[0.75 0.05 0.12 0.05],'Callback','close','String','Close');
set(close_button,'Fontsize',FS); 
