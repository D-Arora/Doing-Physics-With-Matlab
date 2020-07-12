function f = fitFunction(a, x, eqType)

% Evaluate function 
switch eqType
    
    case 1
    % 1: y = a1 * x + a2     linear
    f = a(1) .* x + a(2);
    
    case 2
    % 2: y = a1 * x     linear - proportional
    f = a(1) .* x;    
    
    case 3
    % 3: y = a1 * x^2 + a2 * x + a3     quadratic
    f = a(1) * x.^2 + a(2) .* x + a(3);
    
    case 4
    % 4: y = a1 * x^2     simple quadratic
    f = a(1) .* x.^2;
    
    case 5
    % 5: y = a1 * x^3 + a2 * x^2 + a3 * x + a4       cubic polynomial
    f = a(1) .* x.^3 + a(2) .* x.^2 + a(3) .* x + a(4);
     
    case 6
    % 6: y = a1 * x^a2     power
    f = a(1) .* (x.^a(2));
    
    case 7
    % 7: y = a1 * exp(- a2 * x)     exponential decay
    f = a(1) .* exp(- a(2) .* x);     
    
    case 8
    % 8: y = a1 * (1 - exp(-a2 * x))     "exponetial increase"
    f = a(1) .* (1 - exp(-a(2) .* x)) ; 
    
    case 9
    % 9: y = a1 * x^4 + a2 * x^3 + a3 * x^2 + a4 * x + a5      4th order polynomial
    f = a(1) .* x.^4 + a(2) .* x.^3 + a(3) .* x.^2 + a(4) .* x + a(5);
    
    case 10
   % user defined function
   f = a(1) + a(2) .* exp(-abs((x - a(3))./a(4))) ;
   % f = a(1) + a(2) .* exp(-a(3).* x);
    
end

