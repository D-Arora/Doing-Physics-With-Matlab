% bp005.m

clear
close all
clc



 N = 501;
 
 G_basal = zeros(N,6);
 I_basal = zeros(N,6);
 p = zeros(N,6);

% 1 Basal production of glucose
   a1 = 2.1;
   p(:,1) = linspace(0.8*a1,1.2*a1,N);

% 2 Glucose insulin independent glucose utilization
   a2 = 1.0e-3;
 %  p(:,2) = linspace(0.5e-3, 3e-3, N);
   p(:,2) = linspace(0.8*a2, 1.2*a2, N);
% 3 Insulin sensitivity   
   a3 = 3.0e-3;
 %  p(:,3) = linspace(1.0e-3, 5.0e-3, N);
   p(:,3) = linspace(0.8*a3, 1.2*a3, N);
% 4 Sigmoidal function constant  
   a = 1e4;
%   p(:,4) = linspace(0.5e4, 3e4, N);
   p(:,4) = linspace(0.8*a, 1.2*a, N);
  
% 5 Max rate of insulin secretion from pancreas
   Imax = 0.28;
 %  p(:,5) = linspace(0.1, 0.4, N);
   p(:,5) = linspace(0.8*Imax, 1.2*Imax, N);
   
% 6 Rate at which insulin is cleared (liver mainly)
  tHI = 46;              % half-life  [min]
  kI = log(2)/tHI;
  % p(:,6) = linspace(20, 60, N); 
  p(:,6) = linspace(0.8*tHI, 1.2*tHI, N);
  
% Basal glucose and insulin (roots of cubic polynomial)
  % p1 = a2+a3*Imax/kI; p2 = -a1; p3 = a2*a; p4 = -a1*a;
  % z = roots([p1 p2 p3 p4]);
 
   [Gb, Ib] = basal(a1,a2,a3,a,Imax,kI);
   
   Gb
   Ib

 

 
 for cc = 1: 6
 
     switch cc
       case 1   % Basal production of glucose
       for c = 1 : N
         [G_basal(c,cc), I_basal(c,cc)] = basal(p(c,1),a2,a3,a,Imax,kI);
       end
   
       case 2   % Glucose insulin independent glucose utilization
       for c = 1 : N
         [G_basal(c,cc), I_basal(c,cc)] = basal(a1,p(c,2),a3,a,Imax,kI);
       end
       
       case 3    %  Insulin sensitivity   
       for c = 1 : N
         [G_basal(c,cc), I_basal(c,cc)] = basal(a1,a2, p(c,3),a,Imax,kI);
       end
       
       case 4    %  Sigmoidal function constant
       for c = 1 : N
         [G_basal(c,cc), I_basal(c,cc)] = basal(a1,a2,a3,p(c,4),Imax,kI);
       end 
       
       case 5    %  Max rate of insulin secretion from pancreas
       for c = 1 : N
         [G_basal(c,cc), I_basal(c,cc)] = basal(a1,a2,a3,a,p(c,5),kI);
       end 
       
       case 6    %  Rate at which insulin is cleared (liver mainly)
       for c = 1 : N
         k = log(2)/p(c,6);  
         [G_basal(c,cc), I_basal(c,cc)] = basal(a1,a2,a3,a,Imax,k);
       end 
       
   
     end
 end
 % GRAPHICS  ========================================================== 
 
 figure(1)
  set(gcf,'units','normalized');
  set(gcf,'position',[0.1 0.10 0.35 0.7]);
  set(gcf,'color','w');
  FS = 12;
  
  for cc = 1:6
  
   if cc == 1
       txt = 'a_1';
       tx  =  ' Glucose: basal production';
       x = a1; y = Gb; yR = Ib;
   end
   
   if cc == 2
       txt = 'a_2';
       tx  = 'Glucose: clearance (insulin independent)';
        x = a2; y = Gb; yR = Ib;
   end
   
   if cc == 3
       txt = 'a_3';
       tx  = 'Insulin: sensitivity    ';
        x = a3; y = Gb; yR = Ib;
   end
   
   if cc == 4
       txt = 'a';
       tx  = 'Sigmoidal function constant';
       x = a; y = Gb; yR = Ib;
   end
   
   if cc == 5
       txt = 'I{max}';
       tx  = 'Insulin: max secretion rate from pancreas';
       x = Imax; y = Gb; yR = Ib;
   end
   
   if cc == 6
       txt = 't_H_I';
       tx  = 'Insulin: clearance rate';
        x = tHI; y = Gb; yR = Ib;
   end
     
        
  subplot(3,2,cc)   %  #1
  
    yyaxis left
    xP = p(:,cc); yP = G_basal(:,cc);
    plot(xP,yP,'b','linewidth',2)
    hold on
    hPlot = plot(x,y,'bo');
    set(hPlot,'markerfacecolor','b','markersize',7)
    ylim([60 100])
    grid on
    
    ylabel('G_{basal}')
    set(gca,'fontsize',FS)
    
    yyaxis right
    xP = p(:,cc); yP = I_basal(:,cc);
    plot(xP,yP,'r','linewidth',2)
    hPlot = plot(x,yR,'ro');
    set(hPlot,'markerfacecolor','r','markersize',7)
    ylim([6 12])
      
    xlabel(txt)
    ylabel('I_{basal}')
    set(gca,'fontsize',FS)
    title(tx,'fontweight','normal','fontsize',11)
    
    ax = gca;
    ax.YAxis(1).Color = 'b';
    ax.YAxis(2).Color = 'r';
    
  end

%   subplot(3,2,2)   %  #2
%     yyaxis left
%     xP = aC2; yP = G2_basal;
%     plot(xP,yP,'b','linewidth',2)
%    
%     grid on
%     
%     ylabel('G_{basal}')
%     set(gca,'fontsize',FS)
%     
%     yyaxis right
%     xP = aC2; yP = I2_basal;
%     plot(xP,yP,'r','linewidth',2)
%       
%     xlabel('a_2')
%     ylabel('I_{basal}')
%     set(gca,'fontsize',FS)
%     
%     ax = gca;
%     ax.YAxis(1).Color = 'b';
%     ax.YAxis(2).Color = 'r';
%     
    
     
 
  function  [G_basal, I_basal] = basal(a1,a2,a3,a,Imax,kI)
    p1 = a2+a3.*Imax./kI; p2 = -a1; p3 = a2.*a; p4 = -a1.*a;
    z = roots([p1 p2 p3 p4]);
 
    G_basal = z(1);
    I_basal = (Imax/kI)*G_basal^2/(a+G_basal^2);
  end
  
  