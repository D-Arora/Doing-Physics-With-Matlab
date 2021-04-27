% da_CGLE_01.m



clear 
close all
clc

global k1 k2 h nX

% >>>>> INPUTS
  a = -3;
  b = -2;
  xMax = 100;
  nX = 150;
  tMax = 150;
  nT = 450;

  
% SETUP
t = linspace(0,tMax,nT);
dt = t(2) - t(1);
x = linspace(0,xMax,nX);
dx = x(2) - x(1);
h = dx;

k1 = (1 + 1i*a)/dx^2;
k2 = (1+1i*b);


W = 0.1.*exp(2*1i*x);



% CALCULATIONS


figure (1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.1 0.3 0.3]);
  set(gcf,'color','w');

for n = 1:nT-1
   kx1 = h*tDOT(t(n), W);
   kx2 = h*tDOT(t(n) + h/2, W + kx1/2);
   kx3 = h*tDOT(t(n) + h/2, W + kx2/2);
   kx4 = h*tDOT(t(n) + h,   W + kx3);
      
   W = W + (kx1+2*kx2+2*kx3+kx4)/6;
end
   
    

plot(x,real(W))
hold on
plot(x,imag(W))


% FUNCTIONS

function f = tDOT(t,W)
  global k1 k2 h nX
    
  for c = 2:nX-1
    f = k1*(W(c+1) - 2*W(c) + W(c-1)) + W(c);
    f = f - k2*conj(W(c))*W(c)*W(c);
  end
  
end

