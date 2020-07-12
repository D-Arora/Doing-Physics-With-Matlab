% cemB03.m
% 02 June 2016
% Ian Cooper
% School of Physics, University of Sydney
% http://www.physics.usyd.edu.au/teach_res/mp/mphome.htm

% Magnetic field from a set of parallel wires
%    Biot-Savart Law

% SI units used unless stated otherwise

%-------------------------------------------------------------------------%
%-------Square loop Magnetic Field Calculation..--------------------------%
%-------------------------------------------------------------------------%

clc
close all
clear all

%-------------------------------------------------------------------------%
%--Square loop is in the X-Y plane and Magnetic Field is Evaluated -------%
%-------------at every point in the Y-Z plane(X=0)------------------------%
%-------------------------------------------------------------------------%

Nz=51;  % No. of grids in Z-axis
Ny=51;  % No. of grids in Y-axis
S=12;   % Length of the Square Loop
N=(S+1)*4;   % No of grids in the loop ( X-Y plane)
N4=N/4;
a=S/2;   
I=3;    % current in the loop
u0=1;   % for simplicity, u0 is taken as 1 (permitivity)

%-------These are the coordinates of square loop..........................

Xc(1:N4)=-a:1:a;
Xc(N4+1:2*N4)=a;
Xc(2*N4+1:3*N4)=a:-1:-a;
Xc(3*N4+1:N)=-a;

Yc(1:13)=-a;
Yc(14:26)=-a:1:a;
Yc(27:39)=a;
Yc(40:52)=a:-1:-a;

%-------------------------------------------------------------------------

yp(1:51)=-25:1:25; % Y-coordinates of the plane where we are interested
zp(1:51)=0:1:50;% Z-coordinates of the plane where we are interested

Y(1:Ny,1:Nz)=0; % This array is for 1-d to 2-d conversion of coordinates
Z(1:Ny,1:Nz)=0;

for i=1:Ny
    Y(i,:)=yp(i); % all y-coordinates value in 2-d form
end
for i=1:Nz
    Z(:,i)=zp(i);% all z-coordinates value in 2-d form
end

%-------------------------------------------------------------------------

for A=1:Ny  % for all points in the Y-Z plane, find the magnetic field
for B=1:Nz  % from the loop situated in the X-Y plane
    
%-------------------------------------------------------------------------    
% calculate R-vector from the loop(X-Y plane)to Y-Z plane where we are 
% interested to find the magnetic field and also the dl-vector along the
% loop where current is flowing
% (note that R is the position vector pointing from loop (X-Y plane) to
% the point where we need the magnetic field (in this case Y-Z plane).)
% dl is the current element vector which will make up the square loop------
%-------------------------------------------------------------------------

for i=1:N-1
Rx(i)=-0.5*(Xc(i)+Xc(i+1));
Ry(i)=(Y(A,B)-(0.5*(Yc(i)+Yc(i+1))));
Rz(i)=Z(A,B);
dlx(i)=Xc(i+1)-Xc(i);
dly(i)=Yc(i+1)-Yc(i);
end
Rx(N)=-0.5*(Xc(N)+Xc(1));
Ry(N)=(Y(A,B)-(0.5*(Yc(N)+Yc(1))));
Rz(N)=Z(A,B);
dlx(N)=-Xc(N)+Xc(1);
dly(N)=-Yc(N)+Yc(1);

%--------------------------------------------------------------------------
% for all the elements along loop, determing dl cross R -------------------
% dl cross R is the curl of vector dl and R--------------------------------
% XCross is X-component of the curl of dl and R, similarly we got Y and Z- 
%--------------------------------------------------------------------------

for i=1:N
Xcross(i)=dly(i).*Rz(i);
Ycross(i)=-dlx(i).*Rz(i);
Zcross(i)=(dlx(i).*Ry(i))-(dly(i).*Rx(i));
R(i)=sqrt(Rx(i).^2+Ry(i).^2+Rz(i).^2);
end

%-------------------------------------------------------------------------
% this will be the biot savarts law equation------------------------------
%--------------------------------------------------------------------------

Bx1=(I*u0./(4*pi*(R.^3))).*Xcross;
By1=(I*u0./(4*pi*(R.^3))).*Ycross;
Bz1=(I*u0./(4*pi*(R.^3))).*Zcross;
%--------------------------------------------------------------------------
% now we have  magnetic field from all current elements in the form of an
% array named Bx1,By1,Bz1, now its time to add them up to get total
% magnetic field at the desired point from the loop------------------------
%-------------------------------------------------------------------------

BX(A,B)=0;       % Initialize sum magnetic field to be zero first
BY(A,B)=0;
BZ(A,B)=0;

%--------------------------------------------------------------------------
% here we add all magnetic field from different current elements which 
% make up the loop
%--------------------------------------------------------------------------

for i=1:N   % loop over all current elements along loop    
    BX(A,B)=BX(A,B)+Bx1(i);
    BY(A,B)=BY(A,B)+By1(i);
    BZ(A,B)=BZ(A,B)+Bz1(i);
end

%-------------------------------------------------------------------------

end
end

%------------------------------------------------------------------------
%---BX is a null field and only BY and BZ are color mapped...............
%--------------------------------------------------------------------------

figure(1)
plot(Xc,Yc,'linewidth',3)
axis([-20 20 -20 20])
xlabel('X-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('square loop co-ordinates','fontsize',14)
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(1); 
set(fh, 'color', 'white'); 
grid on

figure(2)
lim1=min(min(BZ));
lim2=max(max(BZ));
steps=(lim2-lim1)/100;
contour(zp,yp,BZ,lim1:steps:lim2)
axis([1 50 -25 25])
xlabel('Z-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('BZ component','fontsize',14)
colorbar('location','eastoutside','fontsize',14);
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(2); 
set(fh, 'color', 'white'); 


figure(3)
lim1=min(min(BY));
lim2=max(max(BY));
steps=(lim2-lim1)/100;
contour(zp,yp,BY,lim1:steps:lim2)
axis([1 50 -25 25])
xlabel('Z-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('BY component','fontsize',14)
colorbar('location','eastoutside','fontsize',14);
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(3); 
set(fh, 'color', 'white'); 


figure(4)
quiver(zp,yp,BZ,BY,2)
axis([1 50 -25 25])
xlabel('Z-axis','fontsize',14)
ylabel('Y-axis','fontsize',14)
title('B-field Vector flow','fontsize',14)
h=gca; 
get(h,'FontSize') 
set(h,'FontSize',14)
h = get(gca, 'ylabel');
fh = figure(4); 
set(fh, 'color', 'white'); 

%--------------------------------------------------------------------------