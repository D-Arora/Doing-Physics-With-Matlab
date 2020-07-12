% op_energy_enclosed.m


% energy transmitted through aperture -------------------------------------
tQmin = 0;                    % angle thetaP  [rad]
tQmax = 2*pi;

rQmin = 0;                    % radius rhoP  [m]
rQmax = aQ;           

tQ1 = linspace(tQmin,tQmax,nQ);      % polar coordinates for P
rQ1 = linspace(rQmin,rQmax,nQ);
[tQ, rQ] = meshgrid(tQ1,rQ1);

WQ = ones(nQ,nQ);

fnQ = rQ .* WQ;

energyQ = simpson2d(fnQ,rQmin,rQmax,tQmin,tQmax)

area = pi*aQ^2

% energy received through observation plane

% Energy enclosed within a circle of prescribed radius -------------------
rho = sqrt(vPxx.^2 + vPyy.^2);   % magnitude of radial coordinate 

tMin = 0;                 % angle theta
tMax = 2*pi;

rhoMin = min(min(rho));   % radius rho
rhoMax = max(max(rho));

fn = rho .* WP;            % integration - polar coordinates  rho*drho* dphi 


for c = 1 : nP
radius = vPx(c);                % prescribed radius
fn(rho > radius) = 0;           % function = 0 outside prescribed radius 
area(c) = simpson2d(fn,rhoMin,rhoMax,tMin,tMax);   % performs integration
fn = rho .* WP;    
end
area * (1/2*pi)^2 * (rP1(2)-rP(1))/3 * (tP1(2)-tP1(1))/3

%area = 100 .* area / max(area);   % normalize area to 100 %


