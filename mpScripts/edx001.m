% edx001.m

% 230325

% https://lpsa.swarthmore.edu/Fourier/Series/DerFS.html

% https://www3.nd.edu/~nancy/Math30650/Matlab/Demos/fourier_series/fourier_series.html

clear
close all
clc


%%   cell 1    Section 8
    clear; close all; %clc 
    syms t n

    fn = t*cos(n*t)
    N = pi;
    L = [0,pi] ;

    F = int(fn,t, L)/N

 
%%   cell 2
syms t  n
f = t*cos(n*t)

F = int(f)


%%  CELL # 3
t = linspace(-pi,pi,999)
f = sqrt(pi^2-t.^2);



F = simpson1d(f,-pi,pi)/(2*pi)

%%  CELL # 3 Square wave Fourier Series
  clear
  close all
  clc

  n = 1:2:200;
  L = length(n);
  N = 501;
  t = linspace(-2*pi,2*pi,N);
  f = zeros(501,L);
  b = zeros(L,1);
  for cn = 1 : L
      b(cn) = 4./(pi.*n(cn));
      fn(:,cn) = b(cn).*sin(n(cn)*t);
  end
  
 F = sum(fn,2);

% GRAPHICS
figure(1)
   plot(t./pi,F,'b','linewidth',1.2)
   grid on
   xlabel('t')
   ylabel('f(t)')

figure(2)
   bar(b)

%%     cell # 3B  Fourier series   Square wave
clear; close all; clc

N = 1000;  Nt = 999;
b = zeros(N,1);
f = zeros(Nt,N);
t = linspace(-2*pi,2*pi,Nt);


flagS = 3;

switch flagS


case 1     % SQUARE WAVE
  for k = 1:2:N
    b(k) = 4/(k*pi);
    f(:,k) = b(k).*sin(k*t);
  end
  F = sum(f,2);

case 2    % SAWTOOTH WAVE   
  for k = 1:N
    b(k) = 2*(-1)^(k+1)/k;
    f(:,k) = b(k).*sin(k*t);
  end
  F = sum(f,2);

case 3  % TRIANGULAR  WAVE
  for k = 1:2:N
    b(k) = 4/(k^2*pi);
    f(:,k) = b(k).*cos(k*t);
  end
  F = pi/2 - sum(f,2);  

end
 % GRAPHICS
figure(1)
   plot(t./pi,F,'b','linewidth',1.2)
   grid on
   xlabel('t')
   ylabel('f(t)')

%figure(2)
 %  bar(b)


 %% CELL #4  square wave
clear; close all; clc

%Create a vector t of 1000 evenly spaced values between -2*pi and 2*pi.
t = linspace(-2*pi,2*pi,1000);
% Create a vector Sq that is equal to the square wave on inputs given by the vector t.
Sq = square(t);

y = 0.8.*sin(t);
y(y>0) = 0.8; y(y<0) = -0.8;

%Plot the square wave with t on the horizontal axis.
plot(t,Sq,'blue')
hold on
plot(t,y,'r')
grid on
%The following set the spacing of ticks and labels on the x and y axes
%It sets the labels and increases the fontsize for readability
xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
yticks([-2 -1 0 1 2])
ylim([-2 2])
xlabel('t'); ylabel('Sq(t)')
set(gca,'fontsize',18)

%% CELL 5  MATRIX CREATION FUNCTIONS
clear; close all; clc

% Create time vector
  t = 0: 0.1: 20;

% Compute magnitude vector
  y = sin(t);

% Compute the noise vector, noiseY
  noiseY = randn(1,201);

% Compute the resulting noisy signal, resSignal
  resSignal = y + noiseY;

% Plot graph
  plot(t,y,'b')
  hold on
  plot(t,resSignal,'r')
  xlabel('time')
  ylabel('Magnitude')
  legend('y','resSignal')
  grid on
  set(gca,'fontsize',12)

%% cell 6    fminsearch
  clear; close all; clc
  f = @(x)  x.^3 - 2.*x + 1;
  [yMin, xMin] = fminsearch(f,0)
  x = linspace(-1,2,501);
  y = f(x);
  plot(x,y)

  %% CELL # 7
  clear; close all; clc
  % Create the function(s) you will need to compute the following coefficients.
  p = @(t) t.*sin(61*t);
  q = @(t) abs(t).*cos(77*t);

% Find the coefficient b61 of the sawtooth wave.
% (Don't forget the correct constant multiple.)
  b61 = (1/pi)*integral(p,-pi,pi);

% Find the coefficient a77 of the triangle wave.
%( Don't forget the correct constant multiple.)
  a77 = (1/pi)*integral(q,-pi,pi);

 b61
 a77


%% CELL # 8     even numbers
   clear; close all; clc
% Initialize the vector evenNumbers to be all zeros
  evenNumbers = zeros(1,10);
% Use a for loop to save the kth *even* number between 2 and 20 
% to the kth entry of the vector evenNumbers

for k = 1:10
    evenNumbers(k) = 2*k;
end
for k = 1:10
  fprintf('k = %2.0f  evenNumber = %2.0f  \n', k,evenNumbers(k))
end



%%
clear; close all; clc
% for k = 1:10
% S(k) = 2*k;
% fprintf('k = %2.0f  evenNumber = %2.0f  \n', k,S(k))
% end
k = 1: 10;
S = 2.*k;
for k = 1:10
   S(k) = k^2
   
end

%% CELL #9   integration
clear; close all; clc


fun = @(x) exp(5*x).*exp(-1i*2.*x)

q = integral(fun,-pi,pi)

% Create the function f(x)=1/(x^3−2x−c) with one parameter, c.

fun = @(x,c) 1./(x.^3-2*x-c);
% Evaluate the integral from x=0 to x=2 at c=5.
q = integral(@(x) fun(x,5),0,2)

fun = @(x,y) 1./( sqrt(x + y) .* (1 + x + y).^2 )

ymax = @(x) 1 - x;
q = integral2(fun,0,1,0,ymax)



%% CELL 10 Fourier series - complex representation
  clear 
  % close all
  clc

% Enter period T >>>>>
  T = 2*pi;
% Enter tims span limits
  tMin = -2*T; tMax = 2*T;
% Enter number of time grid points N  >>>>>
  N = 999;
% Enter max number of frequency components nMax
  nMax = 21;
% Time span  
  t = linspace(tMin,tMax,N);
% Create square wave
  y = sin(t);
  y (y>0) = 1; y(y<0) = -1;
 y = sin(t) +  cos(2*t); 

% Fourier Series 
  a0 = (1/T)*simpson1d(y,0,T);
  k = 2*pi*1i/T;
  a = zeros(nMax,1);
  f = zeros(N,nMax); S = f;
  for n = 1: nMax
      fn = y.*exp(-k*n*t);  
      a(n) = (1/T)*simpson1d(fn,0,T);
      f(:,n) = a(n).*exp(k*n*t) + conj(a(n)).*exp(-k*n*t);
  end
 % Function 
  F = real(a0 + sum(f,2));


% GRAPHICS  ===========================================================
figure(1)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.06 0.08 0.25 0.25]);
   set(gcf,'color','w')    
   plot(t,y,'b','linewidth',3)
   hold on
   plot(t,F,'r')
   grid on
   xticks([-4*pi -3*pi -2*pi -pi 0 pi 2*pi 3*pi 4*pi])
   xticklabels({'-4\pi','-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi','4\pi'})
   yticks([-2 -1 0 1 2])
   ylim([-2 2])
   xlabel('t'); ylabel('Sq(t)')
   set(gca,'fontsize',18)

figure(2)
   set(gcf,'units','normalized');
   set(gcf,'position',[0.06 0.45 0.25 0.25]);
   set(gcf,'color','w')
   bar(a.*conj(a))
   xlabel('n')
   ylabel('|c_n|^2')
   grid on


 %% CELL 11
 close all
 clc
 clear

 [signal,Fs] = audioread('audioGuitar1.wav');
signal = signal(:,1);
L = length(signal);
signal = signal(round(L*0.3):round(L*0.5));
clc
% Copy the procedure in the video to get the single sided spectrum
y = fft(signal);
yMag = abs(y);
N = length(yMag);
%Fs = 3000;
f = 0:(Fs/N):(Fs/2);
yMag = yMag(1:length(f));
% This code plots yMag as a function of f
% You can comment out this code if necessary while building your solution
figure
plot(f,yMag,'-o');
title('Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')
ylabel('|c_n(f)|')
set(gca,'fontsize',18)
xlim([0,100])

