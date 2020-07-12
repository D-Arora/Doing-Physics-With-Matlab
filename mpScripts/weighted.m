% weightedFit.m

% Ian Cooper
% Doing Physics with Matlab

% Fitting an equation to a set of measurements with uncertainties
% The m-script calls the functions
%         chi2test.m     fitfunction.m     partDev.m     

close all


% INPUTS --------------------------------------------------------------
% Plot range
%xmin = 0;  %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%xmax = 250;  %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
xmin = input('Min x value for fit   xmin  =  ');
disp('  ')
xmax = input('Max x value for fit   xmax  =  ');

% Graph Titles
%tm = 'Weighted LS Fit';  %<<<<<<<<<<<<<<<<<<<<<<<
%tx = 'X';  %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%ty = 'Y';  %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
disp('  ')
disp('Enter Plot Title and Axes Labelling') 
disp('  ')
tm = input('Plot Title:    ','s');
disp('  ')
tx = input('X axis label:  ','s');
disp('  ')
ty = input('Y axis label:  ','s');
disp('  ')

% Equation Type 
%eqType = 1;  %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
disp('Select the equation to be fitted to the data')
disp('  ')
disp(' 1: y = a1 * x + a2')                      % linear
disp(' 2: y = a1 * x')                           % linear - proportional
disp(' 3: y = a1 * x^2 + a2 * x + a3')           % quadratic
disp(' 4: y = a1 * x^2')                         % simple quadratic
disp(' 5: y = a1 * x^3 + a2 * x^2 + a3 * x + a4')   % cubic polynomial
disp('  ' );
disp(' 6: y = a1 * x^a2')                        % power
disp('    Zero values for either x or y are not acceptable')
disp('   ');
disp(' 7: y = a1 * exp(- a2 * x)')               % exponential decay or growth
disp('    Zero values for y are not acceptable')
disp('   ');
disp(' 8: y = a1 * (1 - exp(-a2 * x))')          % "growth"
disp(' 9: y = a1 * x^4 + a2 * x^3 + a3 * x^2 + a4 * x +a5')   % 4th order polynomial
disp('  ')
eqType = input('Enter a number from 1 to 9:      ');

% SETUP ---------------------------------------------------------------
% Data must be stored in matrix wData: Col 1 X; Col 2 Y; Col 3 dX; Col 4 dY
n = length(wData(:,1));                % number of data points
m = 1;                                 % number of fitting parameters a1, a2, ...
mMax = 5;                              % max number of coefficients

nx = 200;                              % number of x points for plotting
dx = (xmax - xmin)/nx;                 % increment in x for plotting    
xf = xmin : dx : xmax;                 % x data for plotting fitted fucntion

nIterations = 500;                     % number of iterations for minimizing chi2

if eqType == 1, m = 2; end
if eqType == 2, m = 1; end
if eqType == 3, m = 3; end
if eqType == 4, m = 1; end
if eqType == 5, m = 4; end
if eqType == 6, m = 2; end
if eqType == 7, m = 2; end                
if eqType == 8, m = 2; end
if eqType == 9, m = 5; end

% INITIALIZE MATRICES -------------------------------------------------
wData = sortrows(wData);               % sort matrix by increasing xValues

x = wData(:,1);                        % x measurements
y = wData(:,2);                        % y measurements
dy = wData(:,4);                       % dy uncertainties in y values
dx = wData(:,3);                       % dx uncertainties in x values - used only for error bars = 
w = zeros(1,n);                        % 1/dy
clear a
clear da
xx = ones(n,m);                        % matrix used to get start values for matrix a

u = 0.001;                             % intial value for weighting factor

% Calculation of 1/dy for y data
for cc = 1 : n
    if dy(cc) == 0, w(cc) = 1;
    else w(cc) = 1/dy(cc);
    end
end


% SET STARTING VALUES FOR COEFFICIENTS --------------------------------
switch eqType
    
    case 1     % 1: y = a1 * x + a2     linear
       xx(:,1) = x;
       a = xx\y;
         
    case 2     % 2: y = a1 * x     linear - proportional
      a = y(n)/x(n);
  
    case 3     % 3: y = a1 * x^2 + a2 * x + a3     quadratic
      xx(:,1) = x.*x;
      xx(:,2) = x;
      a = xx\y;
         
    case 4     % 4: y = a1 * x^2     simple quadratic
      xx = x.*x;
      a = xx\y;
   
    case 5     % 5: y = a1 * x^3 + a2 * x^2 + a3 * x + a4       cubic polynomial
      xx(:,1) = x.^3;
      xx(:,2) = x.^2;
      xx(:,3) = x;
      a = xx\y;
       
    case 6     % 6: y = a1 * x^a2     power
      a(2,1) = log10(y(1)/y(n)) / log10(x(1)/x(n));
      a(1,1) = y(1) / x(1)^a(2); 
      
    case 7     % 7: y = a1 * exp(- a2 * x)     exponential
      a(2,1) = log(y(n)/y(1)) / (x(1) - x(n));
      a(1,1) = y(1) / exp(-a(2)*x(1));
     
  case 8     % 8: y = a1 * (1 - exp(-a2 * x))     "growth"
      a(1,1) = y(n); a(2,1) = abs(log(1-y(1)/a(1,1))/x(1));
    %  a(1) = 10;
    %  a(2) = 0.31;
  case 9     % 9: y = a1 * x^4 + a2 * x^3 + a3 * x^2 + a4 * x + a5       4th order polynomial
      xx(:,1) = x.^4;
      xx(:,2) = x.^3;
      xx(:,3) = x.^2;
      xx(:,4) = x;
      a = xx\y;
end

D = zeros(n,1);
f = zeros(n,1);
chi2 = 0;
chi2new = 0;
rchi2 = 0;
pdf = zeros(m,m);
B = zeros(n,1);
CUR = zeros(m,m);
MCUR = zeros(m,m);
MCOV = zeros(m,m);
COV = zeros(m,m);
anew = zeros(m,1);
da = zeros(m,1);
sigma = zeros(m,1);


% Calculations --------------------------------------------------------
for cc = 1 : nIterations

f = fitFunction(a, x, eqType);        % fitted function evaluated at data points 

D = w' .* (y - f);                    % weigthed difference matrix

chi2 = D' * D;                        % chi squared

pdf = partDev(a, x, w, eqType,n,m);   % partail derivatives df/da

B = pdf' * D;                         % B matrix

CUR = pdf' * pdf;                     % Curvature matrix
                                       
for c1 =  1 : m                       % Modified Curvature matrix
for c2 = 1 : m
   MCUR(c1,c2) = CUR(c1,c2) / sqrt(CUR(c1,c1) * CUR(c2,c2));
end   
end
                                       
for c1 =  1 : m
   MCUR(c1,c1) = (1+u)*MCUR(c1,c1);
end   
 
MCOV = inv(MCUR);                      % Modified Covariance matrix
                                       
da = MCOV * B;                         % increments in a coefficients

anew = a + da;                         % new coefficients

f = fitFunction(anew, x, eqType);      % fitted function evaluated at data points 

D = w' .* (y - f);                     % weigthed difference matrix

chi2new = D' * D;                      % chi squred

if chi2new < chi2
    u = u/10; 
    a = anew;
    if abs(chi2-chi2new) < 0.0001, break; end
end

if chi2new > chi2, u = u*10; end

end  %cc  Iterations


ch12 = chi2new;
rchi2 = chi2/(n-m);
u = 0;

pdf= partDev(a, x, w, eqType, n, m);
CUR = pdf' * pdf;                                % Curvature matrix
                                         % Modified Curvature matrix
for c1 =  1 : m                          
for c2 = 1 : m
MCUR(c1,c2) = CUR(c1,c2) / sqrt(CUR(c1,c1) * CUR(c2,c2));
end   
end
                                       
                                        
MCOV = inv(MCUR);                      % Modified Covariance matrix
                                       
for c1 =  1 : m                        % Covriance matrix                          
for c2 = 1 : m
   COV(c1,c2) = MCOV(c1,c2) / sqrt(CUR(c1,c1) * CUR(c2,c2));
end   
end

for c = 1 : m
    sigma(c) = sqrt(COV(c,c));         % Uncertainties in the coefficients
end

prob = chi2test(n-m, chi2);            % chi-squared probability

% Goodness-of-Fit
tgof = '   *** Acceptable Fit   ***';
if rchi2 > 2.5, tgof = '   ??? Fit may not be acceptable ???'; end;
if rchi2 < 0.4, tgof = '   ? Fit may be too good ?'; end;

% Display results -----------------------------------------------------
disp('   ');
wData
disp('   ');

if eqType == 1, disp('1: y = a1 * x + a2     linear'); end;
if eqType == 2, disp('2: y = a1 * x    linear proportional'); end;
if eqType == 3, disp('3: y = a1 * x^2 + a2 * x + a3     quadratic'); end;
if eqType == 4, disp('4: y = a1 * x^2     simple quadratic'); end;
if eqType == 5, disp('5: y = a1 * x^3 + a2 * x^2 + a3 * x + a4       cubic polynomial'); end;
if eqType == 6, disp('6: y = a1 * x^a2     power'); end;
if eqType == 7, disp('7: y = a1 * exp(- a2 * x)     exponential decay'); end;
if eqType == 8, disp('8: y = a1 * (1 - exp(-a2 * x))     "exponetial increase"'); end;   
if eqType == 9, disp('9: y = a1 * x^4 + a2 * x^3 + a3 * x^2 + a4 * x + a5     4th order polynomial'); end;

fprintf('No. measurements  =  %0.0g   \n', n);
fprintf('Degree of freedom  =  %0.0g   \n', (n-m));
fprintf('chi2  =  %0.6g   \n', chi2);
fprintf('Reduced chi2  =  %0.6g   \n', rchi2);
fprintf('Probability =  %0.3g   \n', prob);
disp(tgof);
disp('  ');
disp('Coefficients  a1, a2, ... , am');
disp(a);
disp('Uncertainties in coefficient');
disp(sigma);
disp('  ')


% Graphics ----------------------------------------------------------
yf = fitFunction(a, xf, eqType); 

figure(1)
set(gcf,'PaperType','A4');
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.1 0.1 0.5 0.6])
     
plot(x,y,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5);
hold on 

plot(xf,yf,'b','LineWidth',2);
grid on

title(tm,'FontSize',12);
xlabel(tx,'FontSize',12);
ylabel(ty,'FontSize',12);

yerrP = y + dy;
yerrM = y - dy;
xerrP = x + dx;
xerrM = x - dx;

for c = 1 : n
plot([x(c),x(c)],[yerrM(c), yerrP(c)],'k','LineWidth',2);
plot([xerrM(c),xerrP(c)],[y(c), y(c)],'k','LineWidth',2);
end
set(gca,'XLim',[xmin xmax]);
set(gca,'FontSize',12)










 
   



        
