function U = pot(x);

global Cpot U1 x1
U = U1;

if x > x1
    U = Cpot/x;
end

