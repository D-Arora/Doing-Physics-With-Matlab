clear all
close all
clc


N = 500;
for x = 1: N
    t(x) = x;
    y = sin(0.05*x);
    if y > 0, v(x) = 2;
    elseif y<0 v(x) = 0;
    else v = 1;
    end
end
    figure(1)
    plot(t,v,'linewidth',2)
       