% doppler.m
% Ian Cooper
% School of Physics, University of Sydney

close all
clear all
clc

N = 20;         % number of time steps

disp('   The Figure Window will update every second'); 
disp('  ');
disp('   Enter speed of source from 0 to 350 m/s')
disp('  ')
 
v_s = input('   Speed of source (m/s), v  =  ');

v   = 340;        % velocity of sound
%v_s = 400;       % velocity of source

freq = 1000;

R = zeros(N);    % radius of circles

theta = linspace(0,2*pi,200);   % angle for generating circle

x_max = max([v_s*(N+1) v * (N+1)]);

figure(1)
set(gcf,'Units','normalized');
set(gcf,'Position',[0.1 0.1 0.4 0.6]);
set(gcf,'color','w');

text_1 = 'speed of source   {\itv_S}  = ';
text_2 =  num2str(v_s,'%0.0f');
text_3 = '  m/s     ';
text_4 = 'speed of sound   {\itv} = 340 m/s';
text_p = [text_1 text_2 text_3 text_4];
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

plot(x_p, y_p);

end

pause(1)


if t < N

plot(x,0,'o','Markersize',6,'MarkerFaceColor','w','MarkerEdgeColor','w');

for c = 1 : N
R = v * (t-c);
if R < 0, R = 0; end;

x_p = R * cos(theta) + v_s * (c-1) ;
y_p = R * sin(theta);

plot(x_p, y_p,'w');

end
end
end

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
%pause(d + 0.1);             % waiting for sound end


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
%pause(d + 0.1);             % waiting for sound end

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


