
function scxy = simpsonxy_coeff(num)

%num must be odd
%evaluates the Simpson coefficients
%1 4 2 4 ...2 4 1

%evaluates two dimension Simpson coefficients

sc = 2*ones(num,1);
sc(2:2:num-1) = 4;
sc(1) = 1;
sc(num) = 1;

scx = meshgrid(sc,sc);
scxy = ones(num,num);
scxy(2:2:num-1,:) = scx(2:2:num-1,:)*sc(2);
scxy(3:2:num-2,:) = scx(3:2:num-2,:)*sc(3);
scxy(1,:) = sc';
scxy(num,:) = sc';

