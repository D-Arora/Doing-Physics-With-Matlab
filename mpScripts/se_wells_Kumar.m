% se_wells.m
% Ian Cooper
% School of Physics, University of Sydney
% Potential Wells for solving time independent Schrodinger Equation
% You change the Well parameters in the code for each well type
% Energies values are in eV (electron volts)
% Lengths are in nm (nanometers)
% Other units S.I.
% xMin and xMax define range for x-axis
% x1, x2, ...  well parameters - distance (nm)
% U1, U2, ...  well parameters = potential energy (eV)

clear all
close all
clc

num = 801;             % Number of data points (odd number)

% Constants ------------------------------------------------------------
hbar = 1.055e-34;      % J.s
e = 1.602e-19;         % C
me = 9.109e-31;        % kg  electron mass
mp = 1.67252e-27;      % kg  proton mass
mn = 1.67482e-27;        % kg neutron mass
eps0 = 8.854e-12;      % F/m

m = me;                % Mass of particle


Ese = 1.6e-19;                      % Energy scaling factor  
Lse = 1e-9;                         % Length scaling factor
Cse = -hbar^2/(2*m) / (Lse^2*Ese);  % Schrodinger Eq constant       

% Potential well parameters --------------------------------------------
U = zeros(num,1);
U_matrix = zeros(num-2);

% Potential well types -----------------------------------------------
disp('  ');
disp('  ');
disp('   POTENTIAL WELL TYPES - enter 1 or 2 or 3 ... or 7');
disp('  ');
disp('   1: Square well');
disp(' ');
disp('   2: Stepped well (Asymmetrical well)');
disp(' ');
disp('   3: Double well');
disp(' ');
disp('   4: Sloping well');
disp(' ');
disp('   5: Truncated Parabolic well');
disp(' ');
disp('   6: Morse Potential');
disp(' ');
disp('   7: Parabolic fit to Morse Potential');
disp(' ');
disp('   8: Lattice');
disp(' ');
wellType = input('Specify well type: 1, 2, 3, 4, 5, 6, 7, 8    ');
disp(' ');
disp(' ');

% Potential Wells -----------------------------------------------------
switch wellType
   
% square well **********************************************************
case 1
    xMin = -0.1;         % default = -0.1 nm
    xMax = +0.1;         % default = +0.1 nm
    x1 = 0.05;            % 1/2 well width: default = 0.05 nm
    U1 = -400;            % Depth of well (eV): default = -400 eV
    
    % x1 = 1.5e-6;         % code for deuteron
    % xMin = 0;
    % xMax = +20*x1;
    % U1 = -58.5e6;
    %m = mp*mn/(mp+mn);     % reduced mass - deuteron
    
    % U1 = 0;             % code for infinite square well     
     
   x = linspace(xMin,xMax, num);
   for cn = 1 : num
   if abs(x(cn)) <= x1, U(cn) = U1; end; %if;
   end; %for
   s = sprintf('Potential Well: SQUARE');
   
% step well ************************************************************
% Asymmetrical well: comment the line %if x(cn) > x1/2, U(cn) = 0; end
    
    case 2
    xMin = -0.15;        % default = -0.15 nm
    xMax = +0.15;        % default = +0.15 nm
    x1 = 0.1;            % Total width of well: default = 0.1 nm
    x2 = 0.06;           % Width of LHS well: default = 0.06 nm;
    U1 = -400;           % Depth of LHS well; default = -400 eV;
    U2 = -250;           % Depth of RHS well (eV); default = -250 eV;
    
    x = linspace(xMin,xMax, num); 
    for cn = 1 : num
    if x(cn) >= -x1/2,      U(cn) = U1; end;  %if
    if x(cn) >= -x1/2 + x2, U(cn) = U2; end;  %if
    %if x(cn) > x1/2,        U(cn) = 0;  end; %if
        % comment above line to give an asymmetrical well
    end %for
    s = sprintf('Potential Well: STEPPED');
   
% double well ************************************************************   
case 3   
   U1 = -10;      % Depth of LHS well: default = -440 eV
   U2 = -10;      % Depth of RHS well: default = -400 eV
   U3 = 0;       % Depth of separation section: default = 100 eV 
   x1 = 1.475;      % width of LHS well: default = 0.10 nm
   x2 = 1.475;      % Width of RHS well: default = 0.10 nm
   x3 = 0.05;      % Width of separtion section: default = 0.10 nm
   xEnd = 0.0;    % parameters to define range for x-axis: default = 0.05 nm
   
   
   
   xMin = -(xEnd+x1+x3/2);
   xMax = x3/2+x2+xEnd;
   dx = (xMax-xMin)/(num-1);
   x = xMin : dx : xMax;  
        
   for cn = 1 : num
       if x(cn) >= xMin+ xEnd & x(cn) <= x1 + xMin+ xEnd, U(cn) = U1; end;
       if x(cn) >= x3/2 & x(cn) <= x2+x3/2, U(cn) = U2; end;
       if abs(x(cn)) < x3/2, U(cn) = U3; end
   end
   s = sprintf('Potential Well: DOUBLE');
   
% sloping well ***********************************************************
% Enter energies in eV and distances in nm
    case 4
    % Input parameters
    xMin = -0.1;      % default value = -0.1 nm
    xMax = +0.1;      % default value = + 0.1 nm
    U1 = -1200;        % Depth of LHS well: default = -1200 eV;
    U2 = -200;        % Depth of RHS well: default = -200 eV;
    x1 = 0.05;        % 1/2 width of well: default = 0.05 nm;
    
    x = linspace(xMin,xMax, num);
    intercept = (U1+U2)/2;
    slope = (U2-U1)/(2*x1);
    
   for cn = 1 : num
   if abs(x(cn))<= x1, U(cn) = slope * x(cn) + intercept; end;
   end %for    
   s = sprintf('Potential Well: SLOPING');
   
% parabolic ****************************************************************
case 5   
    xMin = -0.2;                  % default = -0. nm
    xMax = +0.2;                  % default = +0.2 nm
    x1 = 0.2;                     % width default = 0.2 nm;
    U1 = -400;                    % well depth default = -400 eV;
    
   x = linspace(xMin,xMax, num); 
   for cn = 1 : num
   if abs(x(cn))<=x1/2,   U(cn) = -(4*U1/(x1*x1))*x(cn)^2+U1; end;
   end %for    
   s = sprintf('Potential Well: Truncated PARABOLIC');
   
% Morse *******************************************************************
case 6   
    U1 = -1200;          % Depth of well (eV): default = -2000
    x0 = 0.08;
    S = 0.075;
   
    xMin = 0;
    xMax = 0.3;
    
    dx = (xMax-xMin)/(num-1);
    x = xMin : dx : xMax;    
    
    for cn = 1 : num
    U(cn) = U1 * (1 - (exp((x0 - x(cn))/S)-1)^2);
    end
    
    s = sprintf('Potential Well: MORSE potential');
 
   
%Parabolic fit to Morse  incomplete
% case 7   
%    disp('Change Setup values for Morse Potential');
%    for c = 1 : num
%       a1 = (U1-U2)/(R2^2);
%       a2 = (U1-U2)*2*R1/(R2^2);
%       a3 = -U1+(U1-U2)*R1^2/R2^2;
%       U(c) = a1*x(c)^2+a2*x(c)+a3;
%    	if U(c) > 400, U(c) = 400; end;   
%    end %for    
%    s = sprintf('Potential Well: PARABOLIC fit to MORSE potential');

% Lattice
case 8 
    num = 5000;
    wellNum = 12;      % number of wells
    U1 = -350;
    x1 = 0.05;
    x2 = 0.075;
    xEnd = 0.05;
    wellDepth = U1.*ones(wellNum,1);   % depth of wells
    wellWidth = x1.*ones(wellNum,1);
    wellSeparation = x2.*ones(wellNum-1,1);
    wellcenter = zeros(wellNum,1);
    xMin = 0;
    xMax = 2*xEnd + sum(wellSeparation) + 0.5*(wellWidth(1)+wellWidth(wellNum));
    x = linspace(xMin,xMax,num);
    dx = (xMax-xMin)/(num-1);
    U = zeros(num,1);
    wellcenter(1) = xMin+xEnd+wellWidth(1)/2;
    
    for cm = 2: wellNum
        wellcenter(cm) = wellcenter(cm-1) + wellSeparation(cm-1);
    end
    
    for cm = 1 : wellNum
    for cn = 1 : num 
        if abs(x(cn)-wellcenter(cm)) <= wellWidth(cm)/2; U(cn) = wellDepth(cm); end
    end
    end
     
    
    %wellDepth = U1 .* ones(numWells,1);
    %wellWidth = R1 .* ones(numWells,1);
    %defect-----------------------
    %wellDepth(5) = U2;
    %wellWidth(5) = R2;
    % defect----------------------
     %  for cn = 1 : numWells
      %     wellCentre(cn) = R1/2+ R1*(cn-1);
       %    for cm = 1: num
        %   if abs(x(cm)-wellCentre(cn)) < R2/2, U(cm) = wellDepth(cn); end
         %  end        
       %end  
   
s = sprintf('Potential Well: Lattice');

% Colvalent bonding   V shaped potentail wells   
    case 9 
    num = 1000;
    U1 = -1000;
    x1 = 0.01;
    xMin = -0.15;
    xMax = -xMin;
    x = linspace(xMin,xMax,num);
    dx = (xMax-xMin)/(num-1);
    U = zeros(num,1);
    Kc = e/(4*pi*eps0*L0);
    %Uc = 2*Kc / (x1/2+eps);
   
    for cn = 1 : num
      U(cn) =  -Kc * (1/abs(x(cn)-x1/2+eps) + 1/abs(x(cn)+x1/2+eps));
    end      
    U = U + U(end);
    
    for cn = 1 : num
       if U(cn) < U1, U(cn) = U1; end
    end    
    %W = R2/2;
        %w(1) = -R1-2*W; w(2) = w(1) + W; w(3) = w(2)+W;
        %w(4) =  R1; w(5) = w(4) + W; w(6) = w(5)+W;
        %slope(1) = U1/(w(2)-w(1)); b(1) = 0 - slope(1)*w(1);
        %slope(2) = -slope(1); b(2) = U1 - slope(2)*w(2);
        %slope(3) = slope(1); b(3) = 0 - slope(3)*w(4);
        %slope(4) = slope(2); b(4) = U1 - slope(4)*w(5);
        %Ua = 0; Ub = 0;
        %xa = x(c)-R2/2; xb = x(c)+R2/2;
   %for cn = 1 : num
        %if x(cn) > w(1) & x(cn) <= w(2), U(cn) = slope(1)*x(cn)+b(1); end
        %if x(cn) > w(2) & x(cn) <= w(3) , U(cn) = slope(2)*x(cn)+b(2); end
        %if x(cn) > w(4) & x(cn) <= w(5), U(cn) = slope(3)*x(cn)+b(3); end
        %if x(cn) > w(5) & x(cn) <= w(6) , U(cn) = slope(4)*x(cn)+b(4); end
    %    if x(cn) > w(1) & x(cn) <= w(2), U(cn) = U1; end
     %   if x(cn) > w(2) & x(cn) <= w(3) , U(cn) = U1; end
      %  if x(cn) > w(4) & x(cn) <= w(5), U(cn) = U1; end
       % if x(cn) > w(5) & x(cn) <= w(6) , U(cn) = U1; end
        
        
        %if abs(xa)<(R1/2),   Ua = (1-2*abs(xa)/R1) * U1; end;
   %if abs(xb)<(R1/2),   Ub = (1-2*abs(xb)/R1) * U1; end;
   %U(c) = Ua+Ub;
   %end  % for cn  
   s = sprintf('Potential Well: Double: V Shaped');
   
    case 10
       xMin = -0.1;         % default = -0.1 nm
    xMax = +0.1;         % default = +0.1 nm
    x1 = 0.05;            % 1/2 well width: default = 0.05 nm
    U1 = -400;            % Depth of well (eV): default = -400 eV
    
    % x1 = 1.5e-6;         % code for deuteron
    % xMin = 0;
    % xMax = +20*x1;
    % U1 = -58.5e6;
    %m = mp*mn/(mp+mn);     % reduced mass - deuteron
    
    % U1 = 0;             % code for infinite square well     
     
   x = linspace(xMin,xMax, num);
   for cn = 1 : num
   if abs(x(cn)) <= x1, U(cn) = U1; end; %if;
   end; %for
   s = sprintf('Potential Well: SQUARE');
   

   
end;   % switch wellType

% Graphics -------------------------------------------------------------\
figure(1);
set(gcf,'Name','Potential Energy','NumberTitle','off')
plot(x,U,'LineWidth',3);
axis([xMin-eps xMax min(U)-50 max(U)+50]);

title(s);
xlabel('x   (nm)');
ylabel('energy   (eV)');
grid on

% Make potential energy matrix
dx = (x(2)-x(1));
dx2 = dx^2;


for cn =1:(num-2)
    U_matrix(cn,cn) = U(cn+1);
end;


