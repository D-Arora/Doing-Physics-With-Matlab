% mw_doppler.m

% Doppler Effect
% Fig. 1: Animation of wave pattern for a moving source
% Sounds for Doppler Frequencies
% Animation can be saved as annimated gif
%      set f_gif = 1;
%      enter name of file: ag_name =        
% Fig. 2:  Doppler Equation: 
%      moving source / stationary observe
%      moving observer / stationary source

% Ian Cooper
% School of Physics, University of Sydney
% email: ian.cooper@sydney.edu.au
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% 170515


close all;
clear all;
clc;

N = 22;         % number of time steps

disp('   The Figure Window will update every 0.1 s'); 
disp('  ');
disp('   Enter speed of source from 0 to 1000 m/s')
disp('  ')
 
v_s = input('   Speed of source (m/s), v  =  ');

v   = 340;        % velocity of sound
%v_s = 400;       % velocity of source

freq = 1000;      % source frequency

R = zeros(N);    % radius of circles

theta = linspace(0,2*pi,200);   % angle for generating circle

x_max = max([v_s*(N+1) v * (N+1)]);

% In front of moving source
   fFront = (freq * v /(v-v_s));
% In rear of moving source
   fRear = freq * v /(v+v_s);

% ========================================================================
%    ANIMATION
% ========================================================================
% flag: 0 animated gif NOT saved or 1 saved /  file name
   f_gif = 0;
   ag_name = 'ag_doppler_500.gif';
% Delay in seconds before displaying the next image  
   delay = 0.1;  
% Frame counter start
   nt = 1; 

% ========================================================================
%    GRAPHICS
% ========================================================================   
figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.35 0.45]);
set(gcf,'color','w');

text_1 = 'speed of source  {\itv_S}  =  {}';
text_2 =  num2str(v_s,'%8.0f');
text_3 = '  m/s  ';
text_4 = '      speed of sound  {\itv} = 340 m/s';
text1 = strcat(text_1,'  ',text_2,text_3,text_4);

t1 = 'f_{rear} = { }';
t2 = num2str(fRear,'%0.0f');
t3 = '   Hz {            }    ';
t4 = 'f_s = {  }';
t5 = num2str(freq,'%0.0f');
t6 = 'f_{front} = { }';
t7 = num2str(fFront,'%0.0f');
if (v_s >= v), t7 ='-----'; end;
t8 = '  Hz';text2 = strcat(t1,t2,t3,t4,t5,t3,t6,t7,t8);

text_p = [{text1},{text2}];
%text_p = [{text_1 text_2 text_3 text_4},{t1 t2 t3 t2 t4 t5}];

title(text_p,'FontSize',14);

hold on

limits = x_max;
axis([-limits limits -limits limits]);

axis equal
axis off
for t = 1:N

x = v_s * (t-1);    

plot(x,0,'o','Markersize',6,'MarkerFaceColor','r','MarkerEdgeColor','r');

for c = 1 : N
R = v * (t-c);
if R < 0, R = 0; end;

x_p = R * cos(theta) + v_s * (c-1) ;
y_p = R * sin(theta);

plot(x_p, y_p,'b','lineWidth',2);

end

pause(0.1)

   if f_gif > 0 
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
  end

if t < N

plot(x,0,'o','Markersize',6,'MarkerFaceColor','w','MarkerEdgeColor','w');

for c = 1 : N
R = v * (t-c);
if R < 0, R = 0; end;

x_p = R * cos(theta) + v_s * (c-1) ;
y_p = R * sin(theta);

plot(x_p, y_p,'w','lineWidth',2.2);
box on
end
end
   
end

% =======================================================================
%  SOUNDS
% ======================================================================

% source sound
text(-0.6*x_max, -1.1*x_max,'source frequency,  {\itf} = 1000  Hz','FontSize',14)
cf = freq;        
sf = 22050;                 % sample frequency (Hz)
d = 5.0;                    % duration (s)
n = sf * d;                 % number of samples
s = (1:n) / sf;             % sound data preparation
s = sin(2 * pi * cf * s);   % sinusoidal modulation

s = s./max(s);

sound(s, sf);               % sound presentation
pause(d + 0.1);             % waiting for sound end


% in front of moving source
cf = (freq * v /(v-v_s));
if cf > 0
text(1*x_max,0,num2str(cf,'%0.0f'),'FontSize',14)

     
sf = 22050;                 % sample frequency (Hz)
d = 5.0;                    % duration (s)
n = sf * d;                 % number of samples
s = (1:n) / sf;             % sound data preparation
s = sin(2 * pi * cf * s);   % sinusoidal modulation

s = s./max(s);

sound(s, sf);               % sound presentation
pause(d + 0.1);             % waiting for sound end

end
if v_s >= v

text(1*x_max,0,'\it BOOOM !!!','FontSize',14)

% make noise
d = 8;
%rand('state',sum(100 * clock));  % initialize random seed
zzz = rng;
noise = randn(1, n);             % Gausian noise
noise = noise / max(abs(noise)); % -1 to 1 normalization

% play noise
disp('WHITE noise');
sound(noise, sf);                % playing sound
%pause(d + 0.1);                  % waiting for sound end


end

% in rear of moving source
cf = freq * v /(v+v_s);
text(-1.2*x_max,0,num2str(cf,'%0.0f'),'FontSize',14)
     
sf = 22050;                 % sample frequency (Hz)
d = 5.0;                    % duration (s)
n = sf * d;                 % number of samples
s = (1:n) / sf;             % sound data preparation
s = sin(2 * pi * cf * s);   % sinusoidal modulation

s = s./max(s);

sound(s, sf);               % sound presentation
%pause(d + 0.1);             % waiting for sound end


% ======================================================================
% DOPPLER FORMULA
% ======================================================================

%%
figure(2)
set(gca,'fontsize',12);
% Source moving towards a stationary observer
   nF = 5000;
   vs = linspace(0,v-10,nF);
   vo = 0;
   fs = freq;
   fo = fs .* v ./ (v - vs);

   plot(vs,fo,'b','lineWidth',2);
   hold on
   plot([v v],[-1000 max(fo)],'k','lineWidth',2);

% Source moving away a stationary observer
   vs = linspace(0,600,nF);
   fo = fs .* v ./ (v + vs);
  plot(vs,fo,'c','lineWidth',2);
  fo = fs .* v ./ (v + vs);

% Observer moving towards stationary source 
  vo = linspace(0,1000,nF);
  fo = fs .* (v+vo) ./ v;
  plot(vo,fo,'r','lineWidth',2);

% Observer away towards stationary source 
  vo = linspace(0,600,nF);
  fo = fs .* (v-vo) ./ v;
  plot(vo,fo,'m','lineWidth',2); 

  
text(350,4200,'v(sound)','fontsize',14);  
hText = text(40,3200,'souce moving','fontsize',14);
set(hText,'color','b');
  
hText = text(400,800,'souce moving','fontsize',14);
set(hText,'color','c');

hText = text(400,1800,'observer moving','fontsize',14);
set(hText,'color','r');

hText = text(40,0,'observer moving','fontsize',14);
set(hText,'color','m');

xlabel('speed of source or observer');
ylabel('observer freq  f_o  [Hz]');
axis([0 600 -1000 5000])
grid on

