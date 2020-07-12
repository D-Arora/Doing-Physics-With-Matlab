% resonance.m
% 
% Numerically integrate second-order ODE: Damped, driven harmonic oscillator
%
% d^x/dt^2 = -a x - b v + A sin(omega*t)
%
% MATLAB-Monkey.com   9/29/2013


omega = 1;      % resonant frequency = sqrt(k/m)
b = 0.2;        % drag coeficient per unit mass
A = 0.1;        % driving amplitude per unit mass
omega0 = 1.2;   % driving frequency

tBegin = 0;     % time begin
tEnd = 80;      % time end

x0 = 0.2;       % initial position
v0 = 0.8;       % initial velocitie

a = omega^2;    % calculate a coeficient from resonant frequency

% Use Runge-Kutta 45 integrator to solve the ODE
[t,w] = ode45(@derivatives, [tBegin tEnd], [x0 v0]);    
x = w(:,1);     % extract positions from first column of w matrix
v = w(:,2);     % extract velocities from second column of w matrix

plot(t,x);

title('Damped, Driven Harmonic Oscillator');
ylabel('position (m)');
xlabel('time (s)');

set(gcf, 'PaperPosition', [0 0 6 3]);
print ('-dpdf','resonance.pdf');
print ('-dpng','resonance.png','-r100');


    % Function defining derivatives dx/dt and dv/dt
    % uses the parameters a, b, A, omega0 in main program but changeth them not
    function derivs = derivatives(tf,wf,a)
        xf = wf(1);            % wf(1) stores x
        vf = wf(2);            % wf(2) stores v
        dxdt = vf;                                     % set dx/dt = velocity
        dvdt = -a * xf - b * vf + A * sin(omega0*tf);  % set dv/dt = acceleration
        derivs = [dxdt; dvdt];  % return the derivatives
    end

