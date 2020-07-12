% zz03

%% constants
clear
close all
clc

I = 100000000; % current in A
mue0 = 4*pi*10^(-7);
xstep = 0.5;
ystep = 0.5;
zstep = 0.5;
nx = 10;
ny = 10;
nz = 10;
RL = 1;
nL = 30;
count = 0;
%%
xm = 0:xstep:nx;    ym = 0:ystep:ny;    zm = 0:zstep:nz;
[X,Y,Z] = meshgrid(xm,ym,zm);
Bx = zeros(size(X));    By = zeros(size(X));    Bz = zeros(size(X));
% create points (xL,yL,zL) for the curent loop, 
% trough  a parametrization
% phi is the parameter
phi = linspace(0,2*pi,nL);
xL = RL*cos(phi);
yL = RL*sin(phi);
zL = zeros(size(xL));
% store the point's of current loop into the array L 
% move the loop in the origin of the coordinate system, by 
% adding a second constant array to L
L = [xL',yL',zL']+0.5*[nx*ones(size(xL')),ny*ones(size(xL')),nz*ones(size(xL'))];
%% main loop
for z=0:zstep:nz % iteration over z
    
    
    for y=0:ystep:ny % iteration over y
        
        
        for x=0:xstep:nx % iteration over x
            
            
            for i=1:1:nL % iteration over current loop elements
                
                count = count+1;
                 
                
                if i == nL
                    dL = L(end,:)-L(1,:); % i tes Leiterelement
                
                    vec = (L(end,:)+L(1,:))*0.5;
                    
                else
                    dL = L(i+1,:)-L(i,:); % i tes Leiterelement
                
                    vec = (L(i+1,:)+L(i,:))*0.5; % Mittelpunkt des iten Leiterelements
                end
                
               
                
               vecR = [x,y,z]; % Ortsvector an dem B berrechnet wird
               
               vecd = vecR-vec; % verbindungsvektor punkt P zu item Feldverursachenden Leiterelement
               
                
                
                
                dB = mue0*I/(4*pi) * cross(dL',vecd') /(norm(vecd)^2); % beitrag des iten LeiterElements zum B feld am ort vecR
                
                current_point = find(X==x & Y==y & Z==z); % find linear index 
                
                Bx(current_point) = Bx(current_point) + dB(1); % access element of this index and add contribution to B dB
                
                By(current_point) = By(current_point) + dB(2);
                
                Bz(current_point) = Bz(current_point) + dB(3);
                
                
                
                %disp(num2str( count/((nL)*(length(Bx(1,1,:))*length(Bx(1,:,1))*length(Bx(:,1,1)))) ));
            end  % End of current loop elements iteration
            
            
        end       % end of x iteration
        
        
        
        
    end  % end of y iteration
       
    
    
    
end  % end of z iteration
%% post processing
%daspect([1 1 1]);
% axis([0 10 0 10 0 10]);
quiver3(X,Y,Z,Bx,By,Bz,5);
hold on;
%plot3(L(:,1),L(:,2),L(:,3));
gca.PlotBoxAspectRatioMode = 'manual';
% stream line
    t = linspace(0,2*pi,30);      R = 0.5; R2 = 2; R3 = 1.5; R4 = 0.1; R5 = 0.2; R6 = 0.3;
    
    sX = [R *cos(t), R2 * cos(t) , R3 * cos(t), R4 * cos(t) , R5 * cos(t), R6 * cos(t)];
    sY = [ R * sin(t) , R2 * sin(t) , R3 * sin(t)  , R4 * sin(t) , R5 * sin(t) , R6 * sin(t)];
    sZ = [zeros(size(t)),zeros(size(t)) ,zeros(size(t)),zeros(size(t)) ,zeros(size(t)) ,zeros(size(t))];
    
     streamline(X,Y,Z,Bx,By,Bz,sX+nx/2,sY+ny/2,sZ+nx/2);
 % streamline(X,Y,Z,Bx,By,Bz,L(1,1),L(1,2),L(1,3)-1);
 % streamline(X,Y,Z,Bx,By,Bz,X,Y,Z);