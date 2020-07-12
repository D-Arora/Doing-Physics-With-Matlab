% se_measurements.m

% %Plot Potential Energy and eigenvalues -----------------------------
% figure(2)
% set(gcf,'Name','Energy','NumberTitle','off')
% plot(x,U,'lineWidth',2);
% xlabel('position x (m)');
% ylabel('energy U, E_n (eV)');
% title(s);
% hold on
% for cn = 1 : n
% h_plot = plot([xMin xMax], [E(cn) E(cn)]);
% set(h_plot, 'lineWidth',2, 'color', [1 0 0]);
% end

% Plot Wave function and Probability function ----------------------------

qn = input('Enter Quantum Number (1, 2, 3, ...), n  =  ');

figure(3)
set(gcf,'Name','Wave Function','NumberTitle','off')
set(gcf,'Color',[1 1 1]);

subplot(3,1,1)
plot(x,U,'LineWidth',2);
grid on
hold on
h_plot = plot([xMin xMax], [E(qn) E(qn)]);
set(h_plot, 'lineWidth',2, 'color', [1 0 0]);
ylabel('PE,   U (eV)');
tm1 = 'Quantum State,  n  =  ';
tm2 = num2str(qn);
tm = [tm1 tm2];
title(tm);
hold off
subplot(3,1,2)
plot(x,-psi(:,qn),'LineWidth',2);
grid on
ylabel('wave function,   \psi ');

subplot(3,1,3);
plot(x,prob(:,qn),'LineWidth',2);
grid on
xlabel('position x (nm)');
ylabel('prob. density,  \psi^2 (nm^{-1}) ');

% Calculate Expectation Values -------------------------------------------
bra = psi(:,qn)';
Psi = psi(:,qn)';

% probability
ket = Psi;
braket = bra .* ket;
Prob = simpson1d(braket,xMin,xMax);

%position x
ket = x .* Psi;
braket = bra .* ket;
xavg = simpson1d(braket,xMin,xMax);


% Graphical output for <x>
figure(4)
set(gcf,'Name','Expectation Value x','NumberTitle','off')
set(gcf,'Color',[1 1 1]);

subplot(3,1,1)
   plot(x,-Psi,'LineWidth',2);
   grid on
   ylabel('wave function');
tm1 = '<x>  =  ';
tm2 = num2str(xavg);
tm3 = '  nm';
tm = [tm1 tm2 tm3];
title(tm);   
subplot(3,1,2)
   plot(x,ket,'LineWidth',2);
   grid on
   ylabel('ket,   \psi ');
subplot(3,1,3);
   plot(x,braket);
   fill(x,braket,'r')
   grid on
   ylabel('bra ket ');
   hold on
   plot([xavg xavg],[min(braket) max(braket)],'lineWidth',3);
xlabel('position x (nm)');
   hold off

   %position x^2
ket = (x.^2) .* Psi;
braket = bra .* ket;
x2avg = simpson1d(braket,xMin,xMax);

% momentum ip      change length units from nm to m
ket = gradient(Psi,x)* hbar / Lse;
braket = bra .* ket;
ipavg = simpson1d(braket,xMin,xMax);

% momentum ip^2    chnage length unit from nm to m
psi1 = gradient(Psi,x);
ket = -gradient(psi1,x) * hbar^2 / Lse^2;
braket = bra .* ket;
ip2avg = simpson1d(braket,xMin,xMax);

%potential energy U
ket = U' .* Psi;
braket = bra .* ket;
Uavg = simpson1d(braket,xMin,xMax);

%kinetic energy K
psi1 = gradient(Psi,x);
ket = - (hbar^2/(2*me))* gradient(psi1,x) / (Lse^2*Ese);
braket = bra .* ket;
Kavg = simpson1d(braket,xMin,xMax);

%total energy
psi1 = gradient(Psi,x);
ket = - (hbar^2/(2*me))* gradient(psi1,x) / (Lse^2*Ese) + (U.*Psi')';	
braket = bra .* ket;
Eavg = simpson1d(braket,xMin,xMax);


% Calculate uncertainites ------------------------------------------------
deltax = sqrt(x2avg - xavg^2) * Lse;         % length units nm to m
deltaip = sqrt(ip2avg + ipavg^2);            % unit N.s
dxdp = (deltax * deltaip)/hbar;              % N.s.m = J.s   (SI unit)

% Display results of calculations in Command Window
disp('  ');
disp('  ');
fprintf('Quantum number, n  =  %0.0g   \n',qn);
fprintf('Energy, E =  %0.6g  eV \n',E(qn));
fprintf('Total Probability = %0.6g   \n',Prob);
fprintf('<x> = %0.6g   nm\n',xavg);
fprintf('<x^2> = %0.6g   nm^2\n',x2avg);
fprintf('<ip> = %0.6g   N.s\n',ipavg);
fprintf('<ip^2> = %0.6g   N^2.s^2\n',ip2avg);
fprintf('<U> = %0.6g   eV\n',Uavg);
fprintf('<K> = %0.6g   eV\n',Kavg);
fprintf('<E> = %0.6g   eV\n',Eavg);
fprintf('<K> + <U> = %0.6g   eV\n',Kavg+Uavg);
disp('  ')
fprintf('deltax = %0.6g  m \n',deltax);
fprintf('delta| = %0.6g   N.s\n',deltaip);
fprintf('(dx dp)/hbar = %0.6g   \n',dxdp);



