% Example of 'magnetic_field' function.
% Computes the magnetic field of a rectangular current loop
% on the X-Y plane, and computes the loop magnetic moment (m).

% Written by Prof. Yoash Levron, Technion, Israel, 2014.

clc;

meu0 = 4*pi*1e-7; % [H/m]  (Henry / meter)
cur = 1;  % [A] loop current

% resolution in the X-Y plane
dx = 0.001;  % [m]
dy = 0.001;  % [m]

% observation region in the X-Y plane
Xmin = 0.15;  % [m]
Xmax = 0.45;  % [m]
Ymin = 0.6;   % [m]
Ymax = 0.8;   % [m]

% define current loop geometry:
epsd = dx/2;  % shift geometry slightly to avoid infinities 
loop_center_x = 0.3 + epsd;  % [m]
loop_center_y = 0.7 + epsd;  % [m]
edge_len = 0.1;  % [m]

left = loop_center_x-(edge_len/2);
right = loop_center_x+(edge_len/2);
up = loop_center_y+(edge_len/2);
down = loop_center_y-(edge_len/2);
p1 = [left up 0];
p2 = [right up 0];
p3 = [right down 0];
p4 = [left down 0];

FROM = [p1 ; p2 ; p3 ; p4];
TO =   [p2 ; p3 ; p4 ; p1];
CUR =  [cur; cur; cur; cur];

% create observation points matrix
Xvec = Xmin:dx:Xmax;  NX = length(Xvec);
Yvec = Ymin:dy:Ymax;  NY = length(Yvec);
usermem = memory;  MaxPossibleArrayDbl=usermem.MaxPossibleArrayBytes / 8;
if (NX*NY > MaxPossibleArrayDbl / 200)
    'Warning - possibly not enough memory. Program terminated. '
    return
end
[X, Y] = meshgrid(Xvec, Yvec);
R = [X(:), Y(:), 0*X(:)];

% compute the field everywhere
Hmat = magnetic_field( FROM, TO, CUR, R );
HmatZ_vec = Hmat(:,3);  % magnetic field in z direction
HmatZ = reshape(HmatZ_vec,NY,NX);

% display field:
figure(1);
H_Z_dB = 20*log10(abs(HmatZ));
imagesc(Xvec,Yvec, H_Z_dB);
xlabel('X [m]');  ylabel('Y [m]');
title('Magnetic field magnitude in dB');
colorbar;

% MAGNETIC MOMENT
R_far_field = 1e3;  % [m]
H_far_field = magnetic_field( FROM, TO, CUR, [R_far_field 0 0] );
H_far_field_z = H_far_field(1,3);
M = abs( H_far_field_z*4*pi*R_far_field^3 );
disp('The magnetic moment should be equal to the loop-area * loop-current');
disp('magnetic moment in A-m^2 :');
Magnetic_Moment_in_Am2 = M


