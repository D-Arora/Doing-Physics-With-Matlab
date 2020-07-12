function prob = chi2test(dof, chi2Max)

 % Function to calculate how chi2 is distibuted
 %   and give prob of a chi2 value exceeding chi2Max

num = 500;                             % number of calculations for area

flag = 10;                             % flag to end calculation

c = 1;                                 % intialize counter
dx = 0.1;                              % increment in chi2 value

% calculate P(chi2) values --------------------------------------------
while flag > 0.0001
    x(c) = (c-1) * dx;
    y(c) = (x(c)/2)^(dof/2-1)*exp(-x(c)/2)/(2*gamma(dof/2));
    if x(c) > 1, flag = y(c); end ;
    c = c + 1 ;
end

% calculate probabilty of chi2 > chi2Max ------------------------------
xF = linspace(chi2Max, x(c-1),num);
yF = (xF./2).^(dof/2-1).*exp(-xF./2)./(2*gamma(dof/2));

prob = 100 * trapz(xF,yF);

chi2_r = chi2Max / dof;                % reduced chi-squared value

if chi2Max > max(x), prob = 0; end;

% Graphics ------------------------------------------------------------
close
figure(99)
set(gcf,'PaperType','A4');
set(gcf,'Units','Normalized')
set(gcf,'Position',[0.5 0.2 0.4 0.4])

bar(xF, yF,2)
hold on
plot(x,y,'LineWidth',2)

% Labeling graph -----------------------------------------------------
tm = 'chi-squared distribtion';
title(tm,'FontSize',12);
xlabel('\chi ^2','FontSize',12);
ylabel('P(x)','FontSize',12);

tL1 = 'dof  = ';
tL2 = num2str(dof);
tL3 = [tL1 tL2];

tL4 = 'chi2Max  = ';
tL5 = num2str(chi2Max);
tL6 = [tL4 tL5];

tL7 = 'chi2_{reduce}  = ';
tL8 = num2str(chi2_r);
tL9 = [tL7 tL8];

tL10 = 'prob %  = ';
tL11 = num2str(prob, 3);
tL12 = [tL10 tL11];

xPos = max(x); 
yPos = max(y);
text(0.6*xPos, 0.95*yPos,tL3,'FontSize',12);
text(0.6*xPos, 0.88*yPos,tL6,'FontSize',12);
text(0.6*xPos, 0.80*yPos,tL9,'FontSize',12);
text(0.6*xPos, 0.70*yPos,tL12,'FontSize',12);

set(gca,'FontSize',12)

