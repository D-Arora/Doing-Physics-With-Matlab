function [sol_p sol_m] = eq_quadratic(a,b,c)

% solution of the quadratic equation
% a*x^2 + b*x +c = 0

sol_p = (-b + sqrt(b^2-4*a*c))/(2*a);

sol_m = (-b - sqrt(b^2-4*a*c))/(2*a);
