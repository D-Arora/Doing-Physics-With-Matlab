%spP6013.m

clear 
clc
close all

F = [0 15.2	28.6 46.1 59.8 75.1];
I = 0:5;

figure(1)
plot(I,F,'bo')
xlabel('current  [ A ]')
ylabel('force F_B  [ N ]')
grid on
set(gca,'fontsize',14)
