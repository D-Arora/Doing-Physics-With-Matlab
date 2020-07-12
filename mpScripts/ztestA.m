% ztestA.m

%****************************************************************
% Program 2: Calculate first and second derivative numerically
% showing how to write differential operator as a matrix
%****************************************************************
% Parameters for solving problem in the interval 0 < x < L
L = 2*pi; % Interval Length
N = 10; % No. of coordinate points
x = linspace(0,L,N)'; % Coordinate vector
dx = x(2) - x(1); % Coordinate step
% Two-point finite-difference representation of Derivative
D=(diag(ones((N-1),1),1)-diag(ones((N-1),1),-1))/(2*dx);
% Next modify D so that it is consistent with f(0) = f(L) = 0
D(1,1) = 0; D(1,2) = 0; D(2,1) = 0; % So that f(0) = 0
D(N,N-1) = 0; D(N-1,N) = 0; D(N,N) = 0; % So that f(L) = 0
% Three-point finite-difference representation of Laplacian
Lap = (-2*diag(ones(N,1),0) + diag(ones((N-1),1),1) ...
+ diag(ones((N-1),1),-1))/(dx^2);
% Next modify Lap so that it is consistent with f(0) = f(L) = 0
Lap(1,1) = 0; Lap(1,2) = 0; Lap(2,1) = 0; % So that f(0) = 0
Lap(N,N-1) = 0; Lap(N-1,N) = 0; Lap(N,N) = 0;% So that f(L) = 0
% To verify that D*f corresponds to taking the derivative of f
% and Lap*f corresponds to taking a second derviative of f,
% define f = sin(x) or choose your own f
f = sin(x);
% And try the following:
Df = D*f; Lapf = Lap*f;
plot(x,f,'b',x,Df,'r', x,Lapf,'g');
axis([0 5 -1.1 1.1]); % Optimized axis parameters
% To display the matrix D on screen, simply type D and return ...
D % Displays the matrix D in the workspace
Lap % Displays the matrix Lap