% testing_01


clear all
close all
clc

x = 1:20
y = sin(x)

ns = 4;

x2 = x(ns:end)
y2 = y(ns:end-1)

y1 = y(1:end-ns)
x1 = x(1:end-ns)

figure(1)
plot(x1,y1)
hold on
plot(x1,y2)



