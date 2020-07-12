% se_super.m
% Ian Cooper
% School of Physics, University of Sydney
% Animation of superposition of 2 stationary states
% Run se_wells.m to select well and well parameters
% Run se_solve to solve time independent SE for selected well

% Set constants ------------------------------------------------------
nt = 40;                    % number of time steps
Nt = 1;                     % number of cycles for animation

qn(1) = 1; qn(2) = 3;       % states to be summed
ac(1) = 0.5; ac(2) = sqrt(1-ac(1)^2);     % coefficients  
En = abs(E(qn) .* e);       % total energy [J]
wn = En ./ hbar;            % angular frequency of vibration [rad/s]
wnMin = abs(wn(1)-wn(2));
Tn = max(2*pi ./ wnMin);    % period of vibration [s]
dt = Tn / nt;               % time step [s]
t = 0;                      % initial time [s]
pMax = max(psi(:,qn(1)));

Eavg_s = -(ac.^2 * En')/e;  % expectation value for total energy [eV]
disp('  ')
disp('Expectation value for total energy (eV)')
disp(Eavg_s)
figure(12)
set(gcf,'Color',[1 1 1]);
PSI(:,1) = psi(:,qn(1)) .* exp(-i * wn(1) * t);
PSI(:,2) = psi(:,qn(2)) .* exp(-i * wn(2) * t);
PSI_s = ac(1) .* PSI(:,1) + ac(2) .* PSI(:,2); 

plot(x,real(PSI_s),'b','lineWidth',2);

M = getframe(gcf);
[im,map] = rgb2ind(M.cdata,512,'nodither');  %RGB to indexed images
im(1,1,1,nt) = 0;

for cn = 1 : nt*Nt
   PSI(:,1) = psi(:,qn(1)) .* exp(-i * wn(1) * t);
   PSI(:,2) = psi(:,qn(2)) .* exp(-i * wn(2) * t);
   PSI_s = ac(1) .* PSI(:,1) + ac(2) .* PSI(:,2); 

   plot(x,real(PSI_s),'b','lineWidth',2);
   hold on
   plot(x,imag(PSI_s),'Color','r','lineWidth',2);
   plot(x,-16 + conj(PSI_s).*PSI_s,'k','lineWidth',2);
   legend('real \Psi','imaginary \Psi','prob. density');
   title_m = ['n_1  =  ', num2str(qn(1)),'    n_2  =  ',num2str(qn(2))]; 
   title(title_m);
   axis off
   axis([xMin xMax -18 8]);
   pause(0.01)
   
   M = getframe(gcf);
   im(:,:,1,cn) = rgb2ind(M.cdata,map,'nodither');
   
   t = t + dt;
   hold off
end

% position  - calculation of expectation value <x> ---------------------

t = 0;
cn = 0;
%num = 1500;
dt = Tn/(500);        % display three periods
%for cn = 1 : num
while t <= 3*Tn
    cn = cn + 1;   
    time(cn) = t;
   PSI(:,1) = psi(:,qn(1)) .* exp(-i * wn(1) * t);
   PSI(:,2) = psi(:,qn(2)) .* exp(-i * wn(2) * t);
   PSI_s = ac(1) .* PSI(:,1) + ac(2) .* PSI(:,2); 
   bra = x .* PSI_s';
   ket = conj(PSI_s);
   braket = ket .* bra';
   xavg_s(cn) = simpson1d(braket',xMin,xMax);
   t = t + dt;
end
figure(13)
set(gcf,'Color',[1 1 1]);
set(gca,'fontsize',14)
plot(time,real(xavg_s),'lineWidth',2);
xlabel('time t (s)','Fontsize',12);
ylabel('<x>  (nm)','fontsize',12);


%  SAVE ANIMATED GIF ======================================================
% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously
ag_name = 'ag_super13.gif';
delay = 0;
imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);

