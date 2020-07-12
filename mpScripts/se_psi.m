% se_psi.m
% Ian Cooper
% School of Physics, University of Sydney
% Graphical display for solutions of Schrodinger Eq. for a given value of n
% Run the files in the order
    % se_wells.m    se_solve.m    se_psi.m 


disp('   ');
disp('   ');
qn = input('Enter Quantum Number (1, 2, 3, ...), n  =  ');

K = E(n) - U;           % kinetic energy  (eV)
EB = -E(qn);            % binding energy  (eV)

h_figure = figure(99);
clf
set(gcf,'Name','Schrodinger Equation: Bound States');
set(gcf,'NumberTitle','off');
set(gcf,'PaperType','A4');
set(gcf','Units','normalized')
set(gcf,'Position',[0.15 0.05 0.7 0.8]);
set(gcf,'Color',[1 1 1]);

axes('position',[0.1 0.6 0.35 0.32]);
[AX, H1, H2] = plotyy(x,psi(:,qn),x,U);
title_m = [s, '   n =   ', num2str(qn),'   E_n = ', num2str(E(qn))] ;
title(title_m,'Fontsize',12);
xlabel('position x (nm)')
ylabel('wave function   psi');
set(get(AX(2),'Ylabel'),'String','U  (eV)')
set(AX(1),'Xlim',[xMin-eps xMax])
set(AX(1),'YColor',[0 0 1])
set(AX(2),'Xlim',[xMin-eps xMax])
set(AX(2),'Ylim',[min(U)-50 max(K)+50])
set(AX(2),'YColor',[1 0 0])
set(H1,'color',[0 0 1],'LineWidth',3);
set(H2,'color',[1 0 0]);
%set(AX(2),'Ytick',[U1:100:0]);
%set(AX(2),'Ytick',[U1 0 -U1]);
line([xMin xMax],[-EB -EB],'LineWidth',2,'Color',[0 0 0],'Parent',AX(2))
line([xMin xMax],[0 0],'Color',[0 0 1],'Parent',AX(1))

axes('position',[0.1 0.2 0.8 0.32]);
[AX, H1, H2] = plotyy(x,prob(:,qn),x,U);
xlabel('position x (nm)')
ylabel('prob density   psi?');
set(get(AX(2),'Ylabel'),'String','U  (eV)')
set(AX(1),'Xlim',[xMin-eps xMax])
set(AX(1),'YColor',[0 0 1])
set(AX(2),'Xlim',[xMin-eps xMax])
set(AX(2),'Ylim',[min(U)-50 max(U)+50])
set(AX(2),'YColor',[1 0 0])
set(H1,'color',[0 0 1],'LineWidth',3);
set(H2,'color',[1 0 0]);
%set(AX(2),'Ytick',[U1:100:50]);
%set(AX(2),'Ytick',[U1 0 -U1]);
line([xMin xMax],[-EB -EB],'LineWidth',2,'Color',[0 0 0],'Parent',AX(2))

axes('position',[0.57 0.6 0.33 0.32]);
plot(x,U,'r','LineWidth',2);
hold on
plot(x(1:75:num),K(1:75:num),'g+');
plot(x,K,'g');
title('E_n(black)   U(red)    K(green)')  ; 
xlabel('position x (nm)');
ylabel('energies  (eV)');
line([xMin xMax],[-EB -EB],'LineWidth',2,'Color',[0 0 0]);
set(gca,'Xlim',[xMin xMax]);
set(gca,'Ylim',[min(U)-50 max(K)+50]);

% Each point represents the location of the particle after
% a measeurment is made on the system 
axes('position',[0.1 0.05 0.8 0.08]);
hold on
num1 = 10000;
axis off

for c = 1 : num1
xIndex = ceil(1+(num-2)*rand(1,1));
yIndex = rand(1,1);
pIndex = max(prob(:,n))*rand(1,1);
   if pIndex < prob(xIndex,qn);
   plot(x(xIndex),yIndex,'s','MarkerSize',2,'MarkerFaceColor','b','MarkerEdgeColor','none');
   end
end
set(gca, 'Xlim',[xMin xMax]);
set(gca, 'Ylim',[0,1]);

hold off










