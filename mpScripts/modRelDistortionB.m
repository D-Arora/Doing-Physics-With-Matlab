% modRelDistortion.m


% DOING PHYSICS WITH MATLAB: 
%   http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm
% Documentation
%   http://www.physics.usyd.edu.au/teach_res/mp/doc/.htm
% Download Scripts
%   http://www.physics.usyd.edu.au/teach_res/mp/mscripts/

% VISUAL PHYSICS ON LINE
%    http://www.physics.usyd.edu.au/teach_res/hsp/sp/spHome.htm

% Ian Cooper  MatlabVisualPhysics@gmail.com

% 191009   Matlab 2019b


clear
%close all
clc

global xM yM zM c v G tS

% INPUTS ==============================================================
tS = 5;
c = 1;
v = 0.9*c;
G = 1/sqrt(1-v^2/c^2);

xM = 1;
yM = 1;
zM = 1;

% xM = [-1 1 1 -1];
% yM = [ 1 1 3  3];
% zM = [ 1 1 1  1];
% 
% tM = zeros(1,4);
% xS = zeros(1,4);
% k  = zeros(1,4);
% yS = yM;
% zS = zM;

% aQ = 1/G^2;
% 
% for cc = 1:6
%   bQ = -2*(xM(cc)/G + v*tS);
%   cQ = v^2*tS^2 + xM(cc)^2/G^2 + 2*xM(cc)*v*tS/G - v^2*(yM(cc)^2 + zM(cc)^2) /c^2;
%   xS(cc) = (-bQ - sqrt(bQ^2 - 4*aQ*cQ))/(2*aQ);   
% end

% for cc = 1:4
%   aQ = c^2-1;
%   bQ = -2*tS*c^2;
%   cQ = c^2*tS^2 - yM(cc)^2 - zM(cc)^2;
%   k(cc) = (-bQ - sqrt(bQ^2 - 4*aQ*cQ))/(2*aQ);   
% end
% 
%   tM = k./G - v.*xM./c^2;
%   xS = G.*(xM + v.*tM);
  
% figure(1)
% 
% xP = [xS xS(1)]; yP = [yS yS(1)];
% 
% plot(xP,yP)
% %xlim([-20,20])
% hold on
% %plot([xP(1) xP(end)],[yP(1) yP(end)])

% fun = @f; 
% x0 = 2; % initial point
% z = fzero(fun,x0)

tM = linspace(-20,10,100);

fn = f(tM);

fun = @f;
z = fzero(fun,[-10 5])


figure(1)
plot(tM,fn)
grid on

function fn = f(tM)
global xM yM zM c v G tS
fn = sqrt(G^2.*(xM+v.*tM).^ 2+ yM^2  +zM^2)/c + G.*(tM+v*xM/c^2)-tS;
end
% fun = @f;
% z = fzero(fun,[-5 5])
% 
% function  z = f(v,G,xM,yM,zM)
% 
%    z = xM^2+xM*yM+ b;
%    
% end


