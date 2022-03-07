% climateA.m

clear
close all
clc

tic

fun = @temperature;
T0 = [280,250];
T = fsolve(fun,T0);
T;

toc


% *****************************************************************

function F = temperature(T)

% PARAMETERS  *****************************************************
% Stefan-Boltamann constant  [W.m^-2.K^-4]
  s = 5.670374419e-8;
% solar constant [W.m^-2]
  S = 1.366e3;
% emissitivity
  e = 1;
% Transmision and reflection coefficents: S surface ' A atmopshere
%    1 short wavelengths / 2 long wavelengths
  tA1 = 0.43;  tA2 = 0.06;
  aA1 = 0.36;  aA2 = 0.31;
  aS = 0.75;% 0.11;
% Non-radiative interaction between surface and atmosphere [ W.m^-2.K^-1]
  c = 2.5;

%   aA1 = 0;aA2 = 0;
%   tA1 = 1; tA2 = 1;

% constants
  a11 = -tA1*(1 - aS)*S/4;
  a12 = c;
  a13 = e*s*(1-aA2);
  a14 = -e*s;
  a21 = -(1 - aA1 - tA1 + aS*tA1)*S/4;
  a22 = -e*s*(1 - tA2 - aA2);
  a23 = 2*e*s;


F(1) = a11 + a12*T(1) - a12*T(2) + a13*T(1)^4 + a14*T(2)^4;
F(2) = a21 - a12*T(1) + a12*T(2) + a22*T(1)^4 + a23*T(2)^4;

end






