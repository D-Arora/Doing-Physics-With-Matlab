clear all ;
close all ;
clc;
x = 0:1:100; %defines the physical space for the model
y = 0:1:100; %( creates the grid to plot solution values on)
V(101 ,101) = 0; % Set all points to zero .
for i=1:101
    V(i ,1) = 10;
end
% Set one boundary to a certain function ( in Volts )
for i=1:101 V(1 ,i) = 0; end% Set one boundary to a certain function ( in Volts ) .
for i=1:101 V(i ,101) = 0; end% Set one boundary to a certain function ( in Volts ) .
for i=1:101 V(101 ,i) = 0; end% Set one boundary to to a certain function ( in Volts ) .

n=1000; %sets the number of iterations , needs to be quite high to get a %good representation of the solution
for i=1:n ; %for loop using the relaxation method to calculate the potential
for ix = 2:100; % all internal point
for iy = 2:100; % (N ? 2 total internal points for each dimension )
V(iy , ix) = (V(iy-1,ix) + V(iy+1,ix) + V(iy , ix-1) + V(iy , ix+1)) /4;
end
end
end
%V ;
% Print out values for the el ectr ic potential (uncomment i f you really % want to see this , but it ' s 1000 numbers)
[ X , Y ] = meshgrid (x , y) ; % plots a mesh of the solution
subplot (2 ,2 ,1) ; surf (X , Y , V) %plots a surface plot of the solution xlabel ( ' x ' ) ; ylabel ( ' y ' ) ; zlabel ( ' Potential in Volts ' ) ;
subplot (2 ,2 ,2) ; contour (X , Y , V)%shows a ( shaded ) contour plot of the solution [ potential ] = contour (V) ; xlabel ( ' x ' ) ; ylabel ( ' y ' ) ; clabel ( potential ) ;
subplot (2 ,2 ,3) ; contour3 (X , Y , V) %shows a contour plot with colored lines xlabel ( ' x ' ) ; ylabel ( ' y ' ) ; zlabel ( ' Potential in Volts ' ) ;
subplot (2 ,2 ,4) ; meshc(X , Y , V) %shows a mesh plot of the solution xlabel ( ' x ' ) ; ylabel ( ' y ' ) ; zlabel ( ' Potential in Volts ' ) ;
xlabel ( ' x ' ) ; ylabel ( ' y ' ) ; zlabel ( ' Potential in Volts ' ) ;
