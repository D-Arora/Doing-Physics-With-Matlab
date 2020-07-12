function h = wm_sinefn(action)

if (nargin==0)
   action = 'start';
end


switch action

case 'start'

% Create plot figure and axes.
fig = figure('Name','Sine Function','Visible','off','NumberTitle','off');
axes('Units','normalized','Position',[0.10 0.40 0.85 0.55]);

% Call local function to create uicontrols.
ucH = makeuicontrols;

set(fig,'Backingstore','off','UserData',ucH,'Visible','on',...
        'HandleVisibility','callback');
    
% -------------------------------------------------------------------
case 'edittext'

case 'fire'
   ucH = get(gcf,'UserData');
 
   box1 = get(ucH.box1,'String');    % amplitude
   if (isempty(box1))
      box1 = 10;
   else
      box1 = str2double(box1);
   end

   box2 = get(ucH.box2,'String');    % period
   if (isempty(box2))
      box2 = 100;
   else
      box2 = str2double(box2);
  end
  
   box4 = get(ucH.box4,'String');    % initial phase phi
   if (isempty(box4))
      box4 = 0;
   else
      box4 = str2double(box4);
   end


% calculate sine function -----------------------------------------------   
n = 400;

A = box1;    % amplitude
T = box2;    % period
tmin = 0; tmax = 100;   % time
phi = box4;   % initial phase angle

t = linspace(tmin,tmax,n);

y = A .* sin((2*pi) .*t ./T + phi*pi);

% Graphics ------------------------------------------------------------   
plot(t,y,'b','LineWidth',2)
axis([0 tmax -A A])
%title(['Range: ',num2str(range),' ft   Flight time: ',num2str(ttime),' s'])
title('{\ity} = {\itA} sin(2\pi {\itt/T} + \phi)','FontSize',14) 
xlabel('time  {\itt}  (s)','FontSize',14)
ylabel('extension  {\ity}  (m)','FontSize',14)
set(gca,'FontSize',14);
axis([0 100 -10 10]);
grid on

end

% -------------------------------------------------------------------

% MAKEUICONTROLS   Local function to create uicontrols. 
%     function ucH = makeuicontrols(newx,xstr,ystr,stats,degree)
%       from copy
function ucH = makeuicontrols
   pos = [0.05 0.05 0.4 0.6];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   box on
   
fcolor = get(gcf,'Color');
dx = 0.16;
fs = 14;   % font size
%input box location / dimensions

xbo(1) = 0.05; xbo(2) = 0.05; xbo(3) = 0.40; xbo(4) = 0.38;
ybo(1) = 0.20; ybo(2) = 0.05; ybo(3) = 0.20; ybo(4) = 0.15;
xb(1)  = 0.15; xb(2)  = 0.15; xb(3) = 0.15; xb(4)  = 0.15;
yb(1)  = 0.09; yb(2)  = 0.09; yb(3) = 0.09; yb(4)  = 0.09;


% box1------------------------------------------------------------------------
field = [xbo(1) ybo(1) xb(1) yb(1)];

h_uicontrol = uicontrol('Style','text','Units','normalized',...
      'Position',field,'String','Amplitude A: 0 to 10',...
	  'BackgroundColor',fcolor,'ForegroundColor','k');
set(h_uicontrol,'Fontsize',fs);
  
field(1) = field(1) + dx;
  
ucH.box1 = uicontrol('Style','edit','Units','normalized','Position',field,...
         'String','10','BackgroundColor','white','UserData',100,...
         'CallBack','wm_sinefn(''edittext'')');
set(ucH.box1,'Fontsize',fs);   


% box2------------------------------------------------------------------------
field = [xbo(2) ybo(2) xb(2) yb(2)];

h_uicontrol = uicontrol('Style','text','Units','normalized',...
      'Position',field,'String','Period T:       1 to 100',...
	  'BackgroundColor',fcolor,'ForegroundColor','k');
set(h_uicontrol,'Fontsize',fs);  
field(1) = field(1) + dx;  

ucH.box2 = uicontrol('Style','edit','Units','normalized','Position',field,...
         'String','100','BackgroundColor','white','UserData',100,...
         'CallBack','wm_sinefn(''edittext'')');
set(ucH.box2,'Fontsize',fs);     


% box4------------------------------------------------------------------------
field = [xbo(4) ybo(4) xb(4) yb(4)];

h_uicontrol = uicontrol('Style','text','Units','normalized',...
      'Position',field,'String','phase: phi/pi',...
	  'BackgroundColor',fcolor,'ForegroundColor','k');
set(h_uicontrol,'Fontsize',fs);  
field(1) = field(1) + dx;      

ucH.box4 = uicontrol('Style','edit','Units','normalized','Position',field,...
         'String','0','BackgroundColor','white','UserData',0,...
         'CallBack','wm_sinefn(''edittext'')');
set(ucH.box4,'Fontsize',fs);   


% buttons------------------------------------------------------------------------     
     
h_fire_button = uicontrol('Style','Pushbutton','Units','normalized',...
               'Position',[0.75 0.15 0.12 0.05],'Callback','wm_sinefn(''fire'')','String','GO ! !');
set(h_fire_button,'Fontsize',fs);            
               

h_close_button = uicontrol('Style','Pushbutton','Units','normalized',...
               'Position',[0.75 0.05 0.12 0.05],'Callback','close','String','Close');
set(h_close_button,'Fontsize',fs);                  
               
               
