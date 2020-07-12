close all;
clear all;
clc;


%%
n1 = 100;
n2 = 2000;
c1 = 1;
c2 = 5000;

m = (n2-n1) / (c2-c1);

b = n2 - m * c2;

A = zeros(c2,1);
n = zeros(c2,1);

for c = c1 : c2
n(c) = 2*round(0.5*(m * c + b))+1;
end




rMax = 1;
rMin = eps;
r = linspace(rMin, rMax, c2);

dr = r(2)-r(1);


for c = c1:c2;
    f = ones(1,n(c));
    A(c) = r(c) * simpson1d(f,0,2*pi);
end


% figure(1)
% for c = c1 : c2
% t = linspace(0,2*pi,n(c)); 
% x = r(c) .* cos(t);
% y = r(c) .* sin(t);
% 
% plot(x,y,'o');
% axis square
% hold on
% end



Atheory = pi*rMax^2
Acompute = sum(A)* dr

ratio = 100*Atheory/Acompute

%%
close all;
clear all;
clc;

x = 1:9; y = x;
[xx,yy] = meshgrid(x,y);

t = linspace(0,2*pi,9);

x1 = x .* cos(

figure(1)
plot(xx,yy,'o')
axis off
xlabel('xQ'); ylabel('yQ');


%%


