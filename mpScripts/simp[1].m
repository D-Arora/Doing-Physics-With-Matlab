function [sum] = simp(a,b,n)

%  Computes an approximation to the integral from a to b of
%  f(x) dx  using Simpson's rule with n points.

dx = (b-a)/n;

sum = 0;
for i=0:n,              % Loop over nodes.
  x = a + i*dx;
  fx = f(x);            % Evaluate integrand at x.
  if i==0 | i==n, 
    sum = sum + fx/6;   % Add 1/6 * Value of integrand at endpoints.
  else
    sum = sum + fx/3;   % Add 1/3 * Value of integrand at interior points.
  end;

  if i < n,
    xmid = x + .5*dx;
    fxmid = f(xmid);
    sum = sum + 2*fxmid/3;     % Add 2/3 * Value of integrand at midpoints.
  end;
end;

sum = sum*dx;        %  Multiply final result by dx.


%%  Replace this function with the one you want to integrate.  %%

function fx = f(x)
fx = cos(x.^2);
return;
