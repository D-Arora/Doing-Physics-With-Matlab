% cemVectorsB.m

% Ian Cooper
% Email      matlabvisualphysics@gmail.com
% Website    https://d-arora.github.io/Doing-Physics-With-Matlab/
% Documentation   https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/cem001.pdf

% Scalar (dot) product and vector (cross) product of vectors
% Vectros specified in Cartseain components
% Inputs: vectors A   B   C
% Outputs:   magnitude of vectors  Amag   Bmag   Cmag
% Outputs:   dot (scalar) product
% Outputs:   cross (vector) products)
% Angles between vectors

clear all
close all
clc

% INPUTS -----------------------------------------------------------------
     A = [1 0 1]
     B = [0 1 1]
     C = [2 4 5]
     
% OUTPUTS ----------------------------------------------------------------
% magnitudes      
      Amag = norm(A)
      Bmag = norm(B)
      Cmag = norm(C)
% dot products
      AdotA = dot(A,A)
      AdotB = dot(A,B)
      BdotA = dot(B,A)
      BdotC = dot(B,C)
      CdotA = dot(C,A)
      
% cross products 
      AA = cross(A,A)
      AB = cross(A,B)
      BA = cross(B,A)
      BC = cross(B,C)
      CA = cross(C,A)
      
% magnitude of cross products
      AAmag = norm(AA)
      ABmag = norm(AB)
      BAmag = norm(BA)
      BCmag = norm(BC)
      CAmag = norm(CA)
      
 % Angle between vectors [radians] and [degrees]    
     ABangle = asin(norm(AB) /(Amag * Bmag))
     ABangle_deg = rad2deg(ABangle)
     
     BAangle = asin(norm(BA) /(Amag * Bmag))
     BAangle_deg = rad2deg(BAangle)
     
     BCangle = asin(norm(BC) /(Bmag * Cmag))
     BCangle_deg = rad2deg(BCangle)
     
     CAangle = asin(norm(CA) /(Cmag * Amag))
     CAangle_deg = rad2deg(CAangle)

 % Triple Products -------------------------------------------------------    
    AdotBC = dot(A,cross(B,C))
    BdotCA = dot(B,cross(C,A))
    CdotAB = dot(C,cross(A,B))
 
    AdotCB = dot(A,cross(C,B))
    BdotAC = dot(B,cross(A,C))
    CdotBA = dot(C,cross(B,A))
 
    AcrossBC1 = cross(A,cross(B,C))
    AcrossBC2 = B .* dot(A,C) - C .* dot(A,B)
    AcrossBC = cross(A,cross(B,C))
    ABcrossC = cross(cross(A,B),C)

