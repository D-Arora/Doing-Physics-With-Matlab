function [area] = energy_enclosed(fn,x,y)

%xx = meshgrid(x,x);
%yy = meshgrid(y,y);
%[xx yy] = meshgrid(x,y);

rho = sqrt(x.^2+y.^2);

tMin = 0;        % angle
tMax = 2*pi;

rhoMin = min(min(rho));
rhoMax = max(max(rho));




%fn(rho > rhoMax/2) = 0;
%fn(rho > 5) = 0;
fn
rho
fn = rho .* fn

area = simpson2d(fn,rhoMin,rhoMax,tMin,tMax);

