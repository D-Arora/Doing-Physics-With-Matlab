
% Coulomb interaction between neg and pos charges
% Load initial positions and velocities
% Uses Tim initial values for 100 particles

close all
clear
tic

nfile = 'dell/a0506/mat/plasma/coul050924.m';
ndate = '24 sep 05';

% first set of particles are the electrons

% INPUTS --------------------------------------------------------------
np = 100;                       % number of particles: must be even
nt = 1e5;                      % number of time steps
ng = 1e2;                      % spacing between points in plot
dt = 1e-18;                  % time step
Ex = 0;                        % applied - external electric field
rdMin = 1e-11 .*ones(np,np);   % soft potential - min dist between any two particles

n1 = 50;                     % first set of particles
n2 = 50;                     % second set of particles

% Graph Title
tm1 = num2str(np/2);
tm3 = '  electrons   ';
tm4 = '  protons     ';
tm = [tm1 tm3 tm1 tm4];

% Load Data from File
%load data_init_pos_100f;
%load data_init_vel_100f;

% Charges 
e = 1.60218e-19;      % charge on particles: electrons first half of particles
q = -e .* ones(1,np);
q(1+n1:100) = e;
%q(101:130) = 2*e; 
qq = meshgrid(q,q);

% Masses
mP = 1.672623e-27;     % mass of particles: electrons first half of particles
mE = 9.10939e-31;
m = mE * ones(1,np);
%m(1+n1:np) = mP;       % protons
m(1+n1:np) = 47.9 * mP;      % Ti +
% Constant 
eo = 8.85419e-12;
k = 1/(4*pi*eo);
A1 = (k .* q) ./m; 

A2 = A1 .* dt^2;

% SETUP ---------------------------------------------------------------
% Electric Field
Fx = (q .*Ex) ./ m;          %force per unit mass (acceleration)
FxT = Fx .* dt^2;

time = 0 : dt : dt * (nt-1);

% Positions & Distances between particles
x = zeros(nt,np);
y = zeros(nt,np);
z = zeros(nt,np);
rd = zeros(np,np);
rd3 = zeros(np,np);
xd = zeros(np,np);
yd = zeros(np,np);
zd = zeros(np,np);
xx = zeros(np,np);
yy = zeros(np,np);
zz = zeros(np,np);

% Avg distance of particles from origin
rdE = zeros(np,1);
rdP = zeros(np,1);
rEavg = zeros(nt,1);
rPavg = zeros(nt,1);

% Velocities
vx = zeros(nt,np);
vy = zeros(nt,np);
vz = zeros(nt,np);

% Energies
K = zeros(nt,np);
Ktot = zeros(nt,1);
Utot = zeros(nt,1);
Etot = zeros(nt,1);
U =  zeros(np,np);
Upetot = zeros(nt,1) ; Uetot = Upetot ; Uptot = Upetot; 
% Forces - accelerations
SSx = zeros(1,np);
SSy = zeros(1,np);
SSz = zeros(1,np);
Sx = zeros(np,np);
Sy = zeros(np,np);
Sz = zeros(np,np);

% INITIAL CONDITIONS --------------------------------------------------
% data from file
L = 1e-8;                 % length dimension

x0 = zeros(1,np);
y0 = zeros(1,np);
z0 = zeros(1,np);

x0 = L.*(rand(1,np)-0.5) .* ones(1,np);
y0 = L.*(rand(1,np)-0.5) .* ones(1,np);
z0 = L.*(rand(1,np)-0.5) .* ones(1,np);

Te = 3e4;
Tp = 3e4;

scale = 4;
nv = 400;

%flag = 1;        % cold electrons   cold protons
flag = 2;        % hot electrons    cold protons
%flag = 3;        % hot electrons    hot protons

flagM = 1;

% SETUP ---------------------------------------------------------------
ne = n - np;            % number of electrons 

me = 9.10939e-31;
k = 1.38066e-23;
A = -0.5*me/(k*Te);

% cold electrons   cold protons
vx0 = zeros(1,n1);
vy0 = zeros(1,n1);
vz0 = zeros(1,n1);
v = ones(1,ne);

% Maxwellian Distribution 
%    hot electrons   cold protons
if flag == 2
% Maxwellian function 
vProb = sqrt(2*k*Te/me);
vMax = scale* vProb;
vm = linspace(0,vMax,nv);
f = vm.^2 .* exp(A .* vm.^2);
fMax = max(f);
f = f./fMax;

figure(1);
plot(vm,f)

% Velocity magnitudes
for c = 1 : ne
      while flagM == 1
          vR = vMax * rand;
          fC = (vR.^2 .* exp(A .* vR.^2)) ./fMax;
          fR = rand;
            if fR < fC
            v(c) = vR;
            flagM = 0;
            end;
      end
flagM = 1;
end;

% Velocity components
t = (2*pi) .* rand(1,ne) .* ones(1,ne);
p =   (pi) .* (rand(1,ne)-0.5) .* ones(1,ne);

vp = v ./ sqrt(1+sin(p).^2);
vx0(1:ne) = vp .* cos(t);
vy0(1:ne) = vp .* sin(t);
vz0(1:ne) = vp .* sin(p);
end


x(1,:) = 1e-8.*(rand(1,np)-0.5) .* ones(1,np);
y(1,:) = 1e-8.*(rand(1,np)-0.5) .* ones(1,np);
z(1,:) = 1e-8.*(rand(1,np)-0.5) .* ones(1,np);

%x(1,:) = x0; y(1,:) = y0; z(1,:) = z0;
% try smaller volume element
Ls = 4;
x(1,:) = x0./Ls; y(1,:) = y0./Ls; z(1,:) = z0./Ls;

vx(1,:) = vx0; vy(1,:) = vy0; vz(1,:) = vz0;

% START  time step 1 --------------------------------------------------  
xx = meshgrid(x(1,:),x(1,:));
yy = meshgrid(y(1,:),y(1,:));
zz = meshgrid(z(1,:),z(1,:));
xd = xx - xx';
yd = yy - yy';
zd = zz - zz';
rd = sqrt(xd.^2 + yd.^2 + zd.^2) + eye(np,np);
rd = rd + rdMin;
rd3 = rd.^3;

% potential energies
q1 = meshgrid(q,q);
q2 = q1 .* q1';
q3 = ones(np,np)-eye(np,np);
q3 = q2 .* q3;
qf = (q3 ./rd);
Utot(1) = k * sum(sum(qf))/2;

clear q1 q2 

U = zeros(np,np);
for c1 = 1 : n1;
        for c2 = c1 : n1;
          if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
 end;
Uetot(1) = sum(sum(U))/e; 

U = zeros(np,np);
for c1 = n1+1 : np;
        for c2 = c1 : np;
          if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
 end;
Uptot(1) = sum(sum(U))/e; 

U = zeros(np,np);
for c1 = 1 : np;
        for c2 = c1 : np;
          if q(c1)*q(c2) < 0
            if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
        end;
 end;
Upetot(1) = sum(sum(U))/e; 

U = zeros(np,np);

% time step 2 ---------------------------------------------------------
Sx = (qq.*xd) ./rd3;
Sy = (qq.*yd) ./rd3;
Sz = (qq.*zd) ./rd3;

% accelerations
ax = -A1 .* sum(Sx') + Fx;
ay = -A1 .* sum(Sy');
az = -A1 .* sum(Sz');

x(2,:) = x(1,:) + vx(1,:) .* dt + (0.5*dt^2) .* ax;
y(2,:) = y(1,:) + vy(1,:) .* dt + (0.5*dt^2) .* ay;
z(2,:) = z(1,:) + vz(1,:) .* dt + (0.5*dt^2) .* az;

clear x0 y0 z0 vx0 vy0 vz0 ax ay az A1

% TIME STEPS > 2  ------------------------------------------------------
for t = 3 : nt
   xx = meshgrid(x(t-1,:),x(t-1,:));
   yy = meshgrid(y(t-1,:),y(t-1,:));
   zz = meshgrid(z(t-1,:),z(t-1,:));
   xd = xx - xx';
   yd = yy - yy';
   zd = zz - zz';
   rd = sqrt(xd.^2 + yd.^2 + zd.^2) + eye(np,np);
   rd = rd + rdMin.*ones(np,np);
   rd3 = rd.^3;
   Sx = (qq.*xd) ./rd3;
   Sy = (qq.*yd) ./rd3;
   Sz = (qq.*zd) ./rd3;
   SSx = -A2 .* sum(Sx');
   SSy = -A2 .* sum(Sy');
   SSz = -A2 .* sum(Sz');
   
   qf = (q3 ./rd);
   Utot(t-1) = k * sum(sum(qf))/2;
 
U = zeros(np,np);
for c1 = 1 : n1;
        for c2 = c1 : n1;
          if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
 end;
Uetot(t-1) = sum(sum(U))/e; 

U = zeros(np,np);
for c1 = n1+1 : np;
        for c2 = c1 : np;
          if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
 end;
Uptot(t-1) = sum(sum(U))/e; 

U = zeros(np,np);
for c1 = 1 : np;
        for c2 = c1 : np;
          if q(c1)*q(c2) < 0
            if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
        end;
 end;
Upetot(t-1) = sum(sum(U))/e; 

U = zeros(np,np);

%       for c1 = 1 : np;
%       for c2 = c1 : np;
%            if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
%       end; 
%       end;
%       Utot(t-1) = sum(sum(U));
        
   x(t,:) = 2.*x(t-1,:) - x(t-2,:) + SSx + FxT;
         
   y(t,:) = 2.*y(t-1,:) - y(t-2,:) + SSy;
   z(t,:) = 2.*z(t-1,:) - z(t-2,:) + SSz;
end % for t

U = zeros(np,np);
for c1 = 1 : n1;
        for c2 = c1 : n1;
          if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
 end;
Uetot(nt) = sum(sum(U))/e; 

U = zeros(np,np);
for c1 = n1+1 : np;
        for c2 = c1 : np;
          if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
 end;
Uptot(nt) = sum(sum(U))/e; 

U = zeros(np,np);
for c1 = 1 : np;
        for c2 = c1 : np;
          if q(c1)*q(c2) < 0
            if c1 ~= c2, U(c1,c2) = (k*q(c1)*q(c2))/rd(c1,c2); end;
        end ;
        end;
 end;
Upetot(nt) = sum(sum(U))/e; 

U = zeros(np,np);

clear xx yy zz xd yd zd rd rd3 Sx Sy Sz SSx SSy SSz U

% ENERGIES ------------------------------------------------------------
% Kinetic energy 
for t = 2 : nt-1;
    vx(t,:) = (x(t+1,:)- x(t-1,:))./(2*dt);
    vy(t,:) = (y(t+1,:)- y(t-1,:))./(2*dt);
    vz(t,:) = (z(t+1,:)- z(t-1,:))./(2*dt);
    K(t,:) = (0.5.*m) .* (vx(t,:).^2 + vy(t,:).^2 + vz(t,:).^2);
end;
K(nt,:) = K(nt-1,:);
K(1,:) = (0.5.*m) .* (vx(1,:).^2 + vy(1,:).^2 + vz(1,:).^2);
Ktot = sum(K');
Ktote1 = sum(K(1,1:n1)')./e;
Ktotp1 = sum(K(1,1+n1:np)')./e;
Ktotef = sum(K(nt,1:n1)')./e;
Ktotpf = sum(K(nt,1+n1:np)')./e;
% Potential energy and potential energy
Utot(nt) = Utot(nt-1);
Etot = Ktot' + Utot;

Etot = Etot./e;
Utot = Utot./e;
Ktot = Ktot'./e;

% Average distance from Origin ----------------------------------------
for t = 1 : nt;
    rdE  = sqrt(x(t,1:n1).^2 + y(t,1:n1).^2 + z(t,1:n1).^2);
 
    rdP  = sqrt(x(t,1+n1:np).^2 + y(t,1+n1:np).^2 + z(t,1+n1:np).^2);
    rEavg(t) = (1/n1) * sum(rdE);
    rPavg(t) = (1/n2) * sum(rdP);
end;
   r = sqrt(x(nt,1:np).^2 + y(nt,1:np).^2 + z(nt,1:np).^2);

clear rdE rdP m q qq A2 Fx FxT

% GRAPHICS ------------------------------------------------------------
% Explosion
figure(1)
set(gcf,'Papertype','A4');
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,9.5,9.5]);
s = 1e6;
%ng = 100;
hold on
for t = 2 : ng : nt-1;
plot3(s.*x(t,1:np), s.*y(t,1:np), s.*z(t,1:np),'k.')
end;
set(gca,'Fontsize',12);
xlabel('{\itx}  ({\mu}m)','FontSize',12)
ylabel('{\ity}  ({\mu}m)','FontSize',12)
zlabel('{\itz}  ({\mu}m)','FontSize',12)
grid on
set(gca,'Xlim',[-0.1 0.1]);
set(gca,'Ylim',[-0.1 0.1]);
set(gca,'Zlim',[-0.1 0.1]);


%% Explosion
figure(2)
set(gcf,'Papertype','A4');
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,9.5,9.5]);
hold on
grid on
s = 1e9;
%ng = 500;
for t = 2 : ng : nt-1;
plot3(s.*x(t,1:np), s.*y(t,1:np), s.*z(t,1:np),'k.')
%plot3(x(t,n1+1:np), y(t,n1+1:np), z(t,n1+1:np),'ko')
hold on
end;
set(gca,'Fontsize',12);

% %pause(0.1)
% % % %axis([-3 3 -3 3])
% axis equal
% end
xlabel('{\itx}  (nm)','FontSize',12)
ylabel('{\ity}  (nm)','FontSize',12)
zlabel('{\itz}  (nm)','FontSize',12)
%title(tm)
%i_P = [1:np/2];
%i_E = [np/2+1:np];

%plot3(x(1,i_E), y(1,i_E), z(1,i_E),'rp','MarkerFaceColor','r');
%plot3(x(nt,i_E), y(nt,i_E), z(nt,i_E),'gp','MarkerFaceColor','g');
%plot3(x(1,i_P), y(1,i_P), z(1,i_P), 'kd','MarkerFaceColor','k')
%plot3(x(nt,i_P), y(nt,i_P), z(nt,i_P),'md','MarkerFaceColor','m')
%grid off
axis([-20 20 -20 20 -20 20]);
set(gca,'box', 'on')


%%
% energies
figure(5)
set(gcf,'Papertype','A4');
set(gca,'Fontsize',10);
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,9.5,9.5]);

ts = 1e12;

plot(time .* ts,Ktot,'k','LineWidth',1);
xlabel('time  (ps)','FontSize',12)
ylabel('energy  (eV)','FontSize',12)
set(gca,'FontSize',12);
%title(tm)
hold on
plot(time .* ts,Utot,'k','LineWidth',1);
plot(time .* ts,Etot,'k','LineWidth',1);
plot(time .* ts,Uetot,'k','LineWidth',1);
plot(time .* ts,Uptot,'k','LineWidth',1);
plot(time .* ts,Upetot,'k','LineWidth',1);
%legend('Ktot','Utot','Etot','Uetot','Uptot','Upetot');

h_text = text(0.08,1200,'{\itU_p}');
set(h_text,'FontSize',12,'BackGroundcolor','w');

h_text = text(0.02,680,'{\itU_e}');
set(h_text,'FontSize',12,'BackGroundcolor','w');

h_text = text(0.01,360,'{\itK}');
set(h_text,'FontSize',12,'BackGroundcolor','w');

h_text = text(0.02,80,'{\itE}');
set(h_text,'FontSize',12,'BackGroundcolor','w');

h_text = text(0.08,-220,'{\itU}');
set(h_text,'FontSize',12,'BackGroundcolor','w');

h_text = text(0.08,-2000,'{\itU_p_e}');
set(h_text,'FontSize',12,'BackGroundcolor','w');


% DISPLAY -------------------------------------------------------------
disp('   ');
disp(nfile);
disp(ndate);
fprintf('time increment dt (s)  =  %.2e  \n', dt);
fprintf('No. electrons  =  %.0f  \n', n1);
fprintf('No. ions  =  %.0f  \n', n2);
fprintf('simultation time (s)  =  %.2e  \n', dt*nt);
disp('   ');
fprintf('Initial avg electron distance from origin (m) =  %.3e \n', rEavg(1));
fprintf('Initial avg ion distance from origin (m) =  %.3e \n', rPavg(1));
fprintf('Final avg electron distance from origin (m) =  %.3e \n', rEavg(nt));
fprintf('Final avg ion distance from origin (m) =  %.3e \n', rPavg(nt));
disp('   ');
fprintf('Ktot1   (eV) =  %.3f \n', Ktot(1));
fprintf('Ktotf   (eV) =  %.3f \n', Ktot(nt));
fprintf('Ktote1  (eV) =  %.3f \n', Ktote1);
fprintf('Ktotef  (eV) =  %.3f \n', Ktotef);
fprintf('Ktotp1  (eV) =  %.3f \n', Ktotp1);
fprintf('Ktotpf  (eV) =  %.3f \n', Ktotpf);
disp('   ');
fprintf('Utot1   (eV) =  %.3f \n', Utot(1));
fprintf('Utotf   (eV) =  %.3f \n', Utot(nt));
fprintf('Uetot1     (eV) =  %.3f \n', Uetot(1));
fprintf('Uetotf     (eV) =  %.3f \n', Uetot(nt));
fprintf('Uptot1     (eV) =  %.3f \n', Uptot(1));
fprintf('Uptotf     (eV) =  %.3f \n', Uptot(nt));
fprintf('Upetot1    (eV) =  %.3f \n', Upetot(1));
fprintf('Upetotf    (eV) =  %.3f \n', Upetot(nt));
disp('   ');
fprintf('Etot1  (eV) =  %.3f \n', Etot(1));
fprintf('Etotf  (eV) =  %.3f \n', Etot(nt));
disp('   ');
fprintf('execution time (min)  =  %.2f  \n', toc/60);

