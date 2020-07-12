% math_4u_complex


%%  ex2_113.docx

close all
clear all
clc

z1 = 2-sqrt(3) * 1i
z2 = 1+sqrt(3) * 1i

z1c = conj(z1)
z2c = conj(z2)

R1 = abs(z1)
R2 = abs(z2)

t1 = angle(z1)
t2 = angle(z2)

z3 = z1+z2
z4 = z1-z2c
z5 = z1+sqrt(3)*1i*z2

z12 = z1*z2
R12 = abs(z12)
t12 = angle(z12)

z1_2 = z1 / z2
R1_2 = abs(z1_2)
t1_2= angle(z1_2)

z1n = z1^-12
z2n = z2^24

R1n = abs(z1n)
t1n = angle(z1n)

R2n = abs(z2n)
t2n = angle(z2n)

%%
clear all
close all
clc

z = -14+0.53*1i
abs(z)
abs(1+z*1i)

abs(1/z+i)

