% nsAC.m


% COMPLEX GINZBURG-LANDAU EQUATION  CGLE
clear
close all
clc


% INPUTS ==============================================================
 


  h = 1;
  dt = 1e-3;
  NX = 101;
  NT = 7000;
  D1 = 1.0; 
  D2 = 0.5;
  
  C2 = -3.5;
  C1 = 0.556;
 % C1 = (D1 - D2) /(7*(D1+D2));
  

  
% SETUP =============================================================== 
% Grid
  xG = 0:h:h*(NX-1);
  [xx, yy] = meshgrid(xG,xG);
  
% Time
  t = 0:dt:dt*(NT-1);
  
% Constants for 2nd derivates   constants < 1 for stability
  K2 = dt*(1+1i*C2);
  K1 = dt*(1+1i*C1)/h^2;
 
  
% Arrays
  W = zeros(NX,NX);
  WR = zeros(NX,NX); WI= WR;
%  W = 1e-4.*randn(NX,NX);
  
  W(50,45:55) = rand(1,11);
  u = zeros(NT,1);

% W = (1i.*xx+yy).*exp(-0.03.*(xx.^2+yy.^2));  
 
 
figure(1)
   pos = [0.5 0.1 0.30 0.35];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
for k = 1:NT
    
    x = 2:NX-1;
    y = 2:NX-1;
    
    W(x,y) = W(x,y) + dt.*W(x,y) - K2*conj(W(x,y))*W(x,y)*W(x,y);
    W(x,y) = W(x,y) + K1*(W(x+1,y) - 2*W(x,y) + W(x-1,y));
    W(x,y) = W(x,y) + K1*(W(x,y+1) - 2*W(x,y) + W(x,y-1));
    
 %   WR(x,y) = real(W(x,y)); WI = imag(W(x,y));
 %   W(x,y) = WR(x,y) + 1i*WI(x,y);
    
    u(k) = W(50,40);
%   for c1 = 2:NX-1
%       W(1,c1) = W(2,c1);
%       W(c1,1) = W(c1,2);
%       W(NX,c1) = W(NX-1,c1);
%       W(c1,NX) = W(c1,NX-1);
%   end
 
  if rem(k,500) == 0
     pcolor(xG,xG,real(W))
     shading interp
     axis square
    Hcolorbar = colorbar;
  %  set(Hcolorbar,'ylim',[-2.5 2.5])
  %  zlim([-2.5 2.5])
  %  txt = sprintf('t = %2.0f   max %2.1f   min  %2.1f', k, max(max(W)), min(min(W)));
  %  title(txt)
  %  set(gca,'fontsize',12)
     pause(0.0001)
  end
end  


figure(2)
   pos = [0.4 0.1 0.25 0.25];
   set(gcf,'Units','normalized');
   set(gcf,'Position',pos);
   set(gcf,'color','w');
   box on
   
   subplot(2,1,1)
   plot(t, real(u),'linewidth',2)
   subplot(2,1,2)
   plot(u,'linewidth',2)
%    hold on
%    plot(t, UT(2,:),'r')
%    plot(t, UT(3,:),'m','linewidth',2)
%    grid on
%    xlabel('t');
%    ylabel('v_{membrane}');
%    set(gca,'fontsize',12)
% 
%   % FUNCTIONS  ========================================================
%   function u2 = U(N, u1, v1, dt, uF, Iext, e)
%    u2 = zeros(N,N);
%    x = 2:N-1;
%    y = 2:N-1;
%    u2(x,y) = u1(x,y) + dt.*( (u1(x,y) - u1(x,y).^3./3 - v1(x,y))./e + Iext(x,y) );
%    u2(x,y) = u2(x,y) + uF.*(u1(x+1,y) - 2*u1(x,y) + u1(x-1,y));
%    u2(x,y) = u2(x,y) + uF.*(u1(x,y+1) - 2*u1(x,y) + u1(x,y-1));
%   end
% 
% function v2 = V(N, u1, v1, dt, vF, g, d)
%   v2 = zeros(N,N);
%   x = 2:N-1;
%   y = 2:N-1;
%     
%   v2(x,y) = v1(x,y) + dt.*( u1(x,y) - g.*v1(x,y) + d);
%   v2(x,y) = v2(x,y) + vF.*(v1(x+1,y) - 2.*v1(x,y) + v1(x-1,y));
%   v2(x,y) = v2(x,y) + vF.*(v1(x,y+1) - 2*v1(x,y) + v1(x,y-1));
% end