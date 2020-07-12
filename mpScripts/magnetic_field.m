function Hmat = magnetic_field( FROM, TO, CUR, R )
% The function computes the magnetic field H induced by a given conductors
% geometry. The geometry is represented by straight conductors ("current sticks").
% Theory of this numerical technique may be found in
% "Electromagnetic fields and Energy" by Hermann A. Haus, page 322.
%
% Written by Prof. Yoash Levron, Technion, Israel, 2014.
%
% FUNCTION INPUTS
%
% The shape of conductors is represented by "current sticks". For example,
% a square conductor is represented by four sticks.
%
% FROM - an array of vector points indicating where each current
% stick starts. FROM(i,:) is a raw vector (x,y,z), indicating a
% point in 3-D space.  Units in meters [m]
%
% TO - same as FROM. Indicating where each current stick ends.
%
% CUR - Column Vector representing the current of each stick. CUR(i)is
% a scalar. Units in Amperes [A]
%
% R - Observation points. An array of vector points in which the magnetic
% field is to be calculated. R(i,:) is a raw vector (x,y,z), indicating a
% point in 3-D space. Units in meters [m]
%
% FUNCTION OUTPUT
%
% Hmat - the magnetic field H at the observation points. Hmat(i,:) is a
% raw vector (Hx,Hy,Hz), indicating the magnetic field vector in Cartesian
% coordinates. Units are [A/m].

%%%%%%%%%%%%%%%%%%%%% COMPUTATION  %%%%%%%%%%%%%%%

Hmat = 0*R;  % initialize field matrix

% For each stick compute the field at all observation points
for ii = 1:length(CUR), 
    from1 = FROM(ii,:);
    to1 = TO(ii,:);
    cur1 = CUR(ii);
    B = repmat(from1,size(R,1),1) - R;
    C = repmat(to1,size(R,1),1) - R;
    A = C - B;
    
    crs = cross(C,A,2);   % cross product   
    temp1 = (cur1/(4*pi))*( dot(A,C,2)./(sum(C.^2,2).^0.5) - ...
        dot(A,B,2)./(sum(B.^2,2).^0.5) )./sum(crs.^2,2);
    
    Hpp = crs.*repmat(temp1,1,3);
    
    Hmat = Hmat + Hpp;
end