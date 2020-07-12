% se_stationary.m
% Ian Cooper
% School of Physics, University of Sydney
% Animation of a stationary state
% Run se_wells.m to select well and well parameters
% Run se_solve to solve time independent SE for selected well

% Set constants ------------------------------------------------------
nt = 40;               % number of time steps
Nt = 1;                % number of cycles for animation

qn = input('Enter Quantum Number (1, 2, 3, ...), n  =  ');

En = abs(E(qn) * e);                    % total energy [J]
wn = En / hbar;                   % angular frequency of vibration [rad/s]
Tn = 2*pi / wn;                   % period of vibration [s]
dt = Tn / nt;                     % time step [s]
t = 0;                            % initial time [s]
pMax = max(psi(:,qn));

figure(12)
set(gcf,'Color',[1 1 1]);

PSI = psi(:,qn) .* exp(-i * wn * t) ./ pMax;
plot(x,real(PSI),'b','lineWidth',2);
M = getframe(gcf);
[im,map] = rgb2ind(M.cdata,256,'nodither');  %RGB to indexed images
im(1,1,1,nt) = 0;

for cn = 1 : nt
   PSI = psi(:,qn) .* exp(-i * wn * t) ./ pMax;
   plot(x,real(PSI),'b','lineWidth',2);
   hold on
   plot(x,imag(PSI),'Color','r','lineWidth',2);
   plot(x,-2.2 + conj(PSI).*PSI,'k','lineWidth',2);
   legend('real \Psi','imaginary \Psi','prob. density');
   title_m = ['n  =  ', num2str(qn)]; 
   title(title_m);
   axis off
   axis([xMin xMax -3 2]);
   pause(0.1)
   
   M = getframe(gcf);
   im(:,:,1,cn) = rgb2ind(M.cdata,map,'nodither');
   
   t = t + dt;
   hold off
end

%  SAVE ANIMATED GIF ======================================================
% im - images to be saved
% map - color map for images
% ag_name - file name for animated gif
% DelayTime - time delay in seconds between viewing images
% LoopCount - animated gif will continuously
ag_name = 'ag_stationary.gif';
delay = 0;
imwrite(im,map,ag_name,'DelayTime',delay,'LoopCount',inf);

