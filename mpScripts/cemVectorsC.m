% cemVectorsC.m
% 10 march 2016
% Ian Cooper
% School of Physics, University of Sydney
% ../mphome.htm

% ROTATION throug angle theta of the XY axes:  XY --> X'Y'
% Vectors specified in Cartseain components
% Inputs: vector V (Cartesian components in XY system) & rotation angle theta
% Variables: rotation matrix Rmatrix with elements R 
% Outputs:   Cartesian coordinates in frame of reference X'Y'

clear all
close all
clc

% INPUTS -----------------------------------------------------------------
     V = [2 3 0];
     theta = 30;
     
% CALCULATIONS   
   R(1,1) = theta;
   R(1,2) = 90 - theta;
   R(1,3) = 90;
   
   R(2,1) = 90 + theta;
   R(2,2) = theta;
   R(2,3) = 90;
   
   R(3,1) = 90;
   R(3,2) = 90;
   R(3,3) = 0;

   Rmatrix = cosd(R);

% Transformation of vector
   Vdash = (Rmatrix * V')';
 
 Rmatrix
 V
 Vdash
 
 