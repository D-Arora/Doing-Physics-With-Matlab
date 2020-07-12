% testSE2D

% %
% % Program to illustrate solution of 2D Time Dependent Schrodinger equation
% % using leapfrog algorithm
% % based on 'Computational Physics' book by N Giordano and H Nakanishi
% % Section 10.5
% % by Kevin Berwick
% %
%
% % Initialise and set up initial waveform
clear
close all
clc

N=200;
% Set up intial wavepacket;
x_0 = 0.20;   %0.25
y_0 = 0.5;
C=10;
sigma_squared=0.005;
k_0 = 50;     %40

% Discretisation parameters
delta_x=1/200;
delta_t = 1e-5;%=0.00001;
%
% Build a mesh for the values of the probability density function
a= linspace(0, 1, N);
%Use meshgrid to calculate the grid matrices for the x- and y-coordinates,
%using same resolution and scales in x and y directions.
[x,y] = meshgrid(a);
% Create a 2D potential cliff
V=zeros(N,N);
%V(:, 100:200)=-1e3;
% % Create a 2D potential wall
% V=zeros(N,N);
% V(:, 100:200)=1e3;
  V(1:85   ,97:102) = 15000;
  V(115:200,97:102) = 15000;
  V(90:110,97:102) = 15000;

% Calculate psi
psi_stationary=C*exp(-(x-x_0).^2/sigma_squared).*exp(-(y-y_0).^2/sigma_squared);
plane_wave = exp(1i*k_0*x);%+1i*k_0*y);
psi_z=psi_stationary.*plane_wave;
% % % Extract the real and imaginary parts of the wavefunction
%
R_initial=real( psi_z);
I_initial=imag( psi_z);
% % Initialise current real and imaginary parts of psi
%
I_current=I_initial;
R_current=R_initial;
%
% % Initial run of Im(psi) to start off leapfrog process;
%
%
%[I_next] = imag_psi2D(N, I_current, R_current, delta_t, delta_x, V);
%
% % % Do the leapfrog!!
% %
for time_step = 1:1000 %2000
% evaluate R at delta_t, 2*delta_t, 3*delta_t.......
% Time is incremented by t=t+delta_t/2 every call;
[R_next]=real_psi_2D(N, R_current, I_current, delta_t, delta_x, V);
R_current=R_next;
% evaluate I at (3/2)*delta_t, (5/2)*delta_t............
% Time is incremented by t=t+delta_t/2 every call;
[I_next] = imag_psi_2D(N, I_current, R_current, delta_t, delta_x, V);
% calculate psi*psi with R(t) and I(t+delta_t/2) and I(t-delta_t/2)
prob_density=R_current.^2+I_next.*I_current;
I_current=I_next;
% Visualise results with movie. Plot every 10 timesteps for speed
if rem(time_step, 10)== 0
surf(x,y, abs(prob_density.^0.05));    
%pcolor(x,y, abs((prob_density+eps).^0.05));
title('Probability density function');
xlabel('x');
ylabel('y');
zlabel('ps*psi');
axis([0 1 0 1 0 2]);
%view(3)
view(44,55);
axis on;
grid on;
%colormap('copper');
light;
%zlim([0 100])
lighting phong;
camlight('left');
shading interp;
colorbar;
drawnow;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % Calculate the imaginary part of the wavefunction at time
% t=t+delta_t/2,, t + 3*delta_t/2 etc
% % given the value at time t. Vectorise for speed.
function [I_next]= imag_psi_2D(N, I_current, R_current, delta_t, delta_x, V)
I_next= zeros(N,N);
s=delta_t/(2*delta_x^2);
x=2:N-1;
y=2:N-1;
% Calculate the imaginary part of the wavefunction at time t=t+delta_t,
% given the value at time t.
I_next(x,y)=I_current(x,y) +s*(R_current(x+1,y)-2*R_current(x,y)+R_current(x-1,y)+R_current(x,y+1)-2*R_current(x,y)+R_current(x,y-1))...
-delta_t*V(x,y).*R_current(x,y);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % Calculate the real part of the wavefunction at time t=t+delta_t,
% t+2*delta_t etc....
% % given the value at time t. Vectorise for speed.
function [R_next]= real_psi_2D(N, R_current, I_current, delta_t, delta_x, V)
R_next= zeros(N,N);
s=delta_t/(2*delta_x^2);
x=2:N-1;
y=2:N-1;
% Calculate the real part of the wavefunction at time t=t+delta_t,
% given the value at time t.
R_next(x,y)=R_current(x,y) - s*(I_current(x+1,y)-2*I_current(x,y)+I_current(x-1,y)+I_current(x,y+1)-2*I_current(x,y)+I_current(x,y-1))...
+delta_t*V(x,y).*I_current(x,y);
end