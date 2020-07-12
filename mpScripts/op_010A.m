% op_010A.m


% Fabry-Perot Interferometer
% Coefficient of finesse F
% Airy function AF

% Ian Cooper
% School of Physics, University of Sydney
% DOING PHYSICS WITH MATLAB: www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation: www.physics.usyd.edu.au/teach_res/mp/doc/op1010A.htm
% Mscripts: www.physics.usyd.edu.au/teach_res/mp/mscripts
% Matlab 2018b  181201

clear 
close all
clc

% =====================================================================
% Coefficient of finesse F as a function of reflectance R
   R = linspace(0,0.9,500);
   F = F_coeff(R);
   
figure(1)
  FS = 12;
  pos = [0.1 0.1 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');
  
  xP = R; yP = F;
  plot(xP,yP,'linewidth',2)
  xlabel('reflectance  R','fontsize',FS)
  ylabel('finesse coeff  F','fontsize',FS)
  set(gca,'fontsize',FS)
  grid on 
  box on  
   

 % ===================================================================== 
 % Airy function
   delta = linspace(0, 4*pi, 500);
   R = [0.2, 0.4, 0.6, 0.8, 0.95];
   AF = zeros(length(delta),length(R));
   contrast = zeros(1,length(R));
   
   for c1 = 1: length(R)
      AF(:,c1) = AiryF(R(c1), delta);
   end
   
   disp(' R     contrast')
   for c1 = 1: length(R)
      contrast(c1) =  (max(AF(:,c1)) - min(AF(:,c1))) / min(AF(:,c1));
      fprintf('%2.2f   %4.0f  \n',R(c1), contrast(c1))
   end
  
 
  
 
   
  figure(2)
  FS = 12;
  pos = [0.1 0.1 0.25 0.3];
  set(gcf,'Units','normalized');
  set(gcf,'Position',pos);
  set(gcf,'color','w');

  xP = delta ./ pi;
  for c1 = 1: length(R) 
    yP = AF(:,c1);
    plot(xP,yP,'linewidth',2)
    hold on
  end
  xlabel('\delta / \pi','fontsize',FS)
  ylabel('Airy function  A_F','fontsize',FS)
  set(gca,'fontsize',FS)
  grid on 
  box on  
  legend('0.2','0.4','0.6','0.8','0.95','location','northoutside', ...
       'orientation','horizontal')

% FUNCTION SECTION ====================================================  
  
function F = F_coeff(R)
    F = 4*R ./ (1 - R).^2;
end

function AF = AiryF(R, delta)
   F = F_coeff(R);
   AF = 1 ./ (1 + F.* sin(delta/2).^2);
end

