function integral = simpson1d(f,a,b)

% [1D] integration - Simpson's 1/3 rule
%      f function    a = lower bound    b = upper bound
%      Must have odd number of data points
%      Simpson's coefficients   1 4 2 4 ... 2 4 1

numS = length(f);               % number of data points

if mod(numS,2) == 1
sc = 2*ones(numS,1);
sc(2:2:numS-1) = 4;
sc(1) = 1; sc(numS) = 1;

h = (b-a)/(numS-1);

integral = (h/3) * f * sc;

else
    
integral = 'Length of function must be an ODD number' 
end


