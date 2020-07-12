function ydot = funct(x,y,flag,E);

global Cse
%k = pi;

%ydot = [y(2); -k^2*(y(1)+ pot(x))];

ydot = [y(2); Cse*(E-pot(x))*y(1)];

%ydot = [y(2); Cse*(E*y(1))];