% se_orthonormal.m

% Input two quantum numbers for two states
%    to test normalization or that two different 
%    eigenvectors are orthogonal to each other orthogonal


% Inputs -----------------------------------------------------
disp('   ');
qn1 = input('Enter Quantum Number (1, 2, 3, ...), n1  =  ');
disp('   ');
qn2 = input('Enter Quantum Number (1, 2, 3, ...), n2  =  ');
disp('   ');

% Evaluate integral -------------------------------------------
bra = psi(:,qn1)';
ket = psi(:,qn2)';
braket = bra .* ket;
integral = round(simpson1d(braket,xMin,xMax));

% Graphical output for evaluation of integral ----------------
figure(4)
set(gcf,'Name','Orthonormal','NumberTitle','off')
set(gcf,'Color',[1 1 1]);

subplot(3,1,1)
   plot(x,bra,'LineWidth',3);
   grid on
   ylabel('\psi_1');
tm1 = 'integral   =  ';
tm2 = num2str(integral,2);
tm = [tm1 tm2];
title(tm);   
subplot(3,1,2)
   plot(x,ket,'LineWidth',3);
   grid on
   ylabel('\psi_2');
subplot(3,1,3);
   plot(x,braket);
   fill(x,braket,'r')
   grid on
   ylabel('bra ket ');
   hold on
   xlabel('position x (nm)');
   hold off

% Display results of calculations in Command Window --------------------
fprintf('Integral =  %0.1f  \n',integral);
disp('   ');
if integral < 0.5
    disp('Wavefunctions are orthogonal');
else
    disp('Wavefunction is normalized');
end




