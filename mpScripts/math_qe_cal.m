% math_qe_cal.m



pInitial = uicontrol(gcf,'Style','text','Position',[350 300 350 110], ...
    'String','   ','BackgroundColor',[0.95 0.9 0.9]);

%pInitial = uicontrol(gcf,'Style','text','Position',[100 150 350 110], ...
%    'String','   ','BackgroundColor',[0.95 0.9 0.9]);

% Impossible equation ----------------------------------------------------
if a ==0 && b == 0;
   p(1) = uicontrol(gcf,'Style','text','Position',[350 380 350 30],...
          'String','This is impossible equation','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]);
   
   p(14) = uicontrol(gcf,'Style','text','Position',[100 200 200 30],...
          'String',' Turning point (x) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
   p(15) = uicontrol(gcf,'Style','text','Position',[300 200 140 30],...
          'String','---','FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9], 'ForegroundColor',[0 0 1]); 
      
   p(16) = uicontrol(gcf,'Style','text','Position',[100 150 200 30],...
          'String','Turning point (y) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
   p(17) = uicontrol(gcf,'Style','text','Position',[300 150 140 30],...
          'String','---','FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9],'ForegroundColor',[0 0 1]);  
flag = 0;
end 

% Linear equation --------------------------------------------------------
if a == 0 && b ~= 0;
    xL = - c / b;
    
    p(2) = uicontrol(gcf,'Style','text','Position',[350 380 350 30],...
          'String','This is a linear equation','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]);
      
    p(3) = uicontrol(gcf,'Style','text','Position',[350 340 170 30],...
          'String',' y = 0  at  x =  ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]);  
      
    p(4) = uicontrol(gcf,'Style','text','Position',[495 340 165 30],...
          'String',xL,'FontSize',14,'HorizontalAlignment','Left',...
          'BackgroundColor',[0.95 0.9 0.9], 'ForegroundColor',[1 0 0]);  
      
   xPmin = xL - 5*xL;
   xPmax = xL + 5*xL;
   xP = linspace(xPmin,xPmax,500);
   yP = b * xP + c;
   
       p(14) = uicontrol(gcf,'Style','text','Position',[100 200 200 30],...
          'String',' Turning point (x) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(15) = uicontrol(gcf,'Style','text','Position',[300 200 140 30],...
          'String','---','FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9], 'ForegroundColor',[0 0 1]); 
      
        p(16) = uicontrol(gcf,'Style','text','Position',[100 150 200 30],...
          'String','Turning point (y) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(17) = uicontrol(gcf,'Style','text','Position',[300 150 140 30],...
          'String','---','FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9],'ForegroundColor',[0 0 1]);  

flag = 1;
end

% Quadrtatic case ---------------------------------------------------------

if a ~= 0;
    
Del = b^2 - 4 * a * c;

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     if Del < 0
         p(5) = uicontrol(gcf,'Style','text','Position',[350 380 350 30],...
          'String','No real roots','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]);
      
       p(14) = uicontrol(gcf,'Style','text','Position',[100 200 200 30],...
          'String',' Turning point (x) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(15) = uicontrol(gcf,'Style','text','Position',[300 200 140 30],...
          'String','---','FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9], 'ForegroundColor',[0 0 1]); 
      
        p(16) = uicontrol(gcf,'Style','text','Position',[100 150 200 30],...
          'String','Turning point (y) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(17) = uicontrol(gcf,'Style','text','Position',[300 150 140 30],...
          'String','---','FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9],'ForegroundColor',[0 0 1]);   
     flag = 0;
     end
 % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     if Del == 0;
        x = - b / (2*a);
        
        p(6) = uicontrol(gcf,'Style','text','Position',[350 380 350 30],...
          'String','There is only one real root','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(7) = uicontrol(gcf,'Style','text','Position',[350 340 50 30],...
          'String',' x = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(8) = uicontrol(gcf,'Style','text','Position',[400 340 165 30],...
          'String',x,'FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9], ...
          'ForegroundColor',[1 0 0]); 
     
     x_turn = - b / (2*a);
     xPmin = x_turn - 10*x_turn;
     xPmax = x_turn + 10*x_turn;
          
     if b == 0;
        xPmin = -10; xPmax = 10; 
     end
     x_m = x_turn;
     y_m = a * x_m^2 + b * x_m + c;
     
     xP = linspace(xPmin,xPmax,500);
     yP = a * xP.^2 + b .* xP + c;
     
     flag = 1;
     
        p(14) = uicontrol(gcf,'Style','text','Position',[100 200 200 30],...
          'String',' Turning point (x) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(15) = uicontrol(gcf,'Style','text','Position',[300 200 140 30],...
          'String',x_m,'FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9], 'ForegroundColor',[0 0 1]); 
      
        p(16) = uicontrol(gcf,'Style','text','Position',[100 150 200 30],...
          'String','Turning point (y) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(17) = uicontrol(gcf,'Style','text','Position',[300 150 140 30],...
          'String',y_m,'FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9],'ForegroundColor',[0 0 1]);   
        
     
     end
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^     
     if Del > 0;
        x1 = (-b + sqrt(Del)) / (2*a);
        x2 = (-b - sqrt(Del)) / (2*a);
        
        p(9) = uicontrol(gcf,'Style','text','Position',[350 380 350 30],...
          'String','There are two real root','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
       
        p(10) = uicontrol(gcf,'Style','text','Position',[350 340 50 30],...
          'String',' x1 =    ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(11) = uicontrol(gcf,'Style','text','Position',[405 340 165 30],...
          'String',x1,'FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9], 'ForegroundColor',[1 0 0]); 
      
        p(12) = uicontrol(gcf,'Style','text','Position',[350 300 50 30],...
          'String',' x2 =    ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(13) = uicontrol(gcf,'Style','text','Position',[405 300 165 30],...
          'String',x2','FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9],'ForegroundColor',[1 0 0]); 
          
        x_m = -b / (2*a);
        y_m = a * x_m^2 + b * x_m + c;
        
        p(14) = uicontrol(gcf,'Style','text','Position',[100 200 200 30],...
          'String',' Turning point (x) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(15) = uicontrol(gcf,'Style','text','Position',[300 200 140 30],...
          'String',x_m,'FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9], 'ForegroundColor',[0 0 1]); 
      
        p(16) = uicontrol(gcf,'Style','text','Position',[100 150 200 30],...
          'String','Turning point (y) = ','FontSize',14,...
          'HorizontalAlignment','Left', 'BackgroundColor',[0.95 0.9 0.9]); 
      
        p(17) = uicontrol(gcf,'Style','text','Position',[300 150 140 30],...
          'String',y_m,'FontSize',14,'HorizontalAlignment','Left', ...
          'BackgroundColor',[0.95 0.9 0.9],'ForegroundColor',[0 0 1]);   
        
     x_turn = - b / (2*a);
     xPmin = x_turn - 10*x_turn;
     xPmax = x_turn + 10*x_turn;
     
     if b == 0;
        xPmin = -10; xPmax = 10; 
     end
     
     xP = linspace(xPmin,xPmax,500);
     yP = a * xP.^2 + b .* xP + c;
     
     flag = 2;
     end

end  
 
     if flag == 0 ;
     xP = zeros(1,500); yP = zeros(1,500);
     end
         
     plot3 = subplot('Position',[0.6 0.1 0.35 0.35]);
     plot(xP,yP,'k','lineWidth',2);
     axis on; grid on;
    
     if flag == 2;
         hold on
         hp1 = plot(x_m,y_m,'o');
         set(hp1,'MarkerSize',6,'MarkerEdgeColor',[0 0 1], 'MarkerFacecolor', [0 0 1]);
         hp1 = plot(x1,0,'o');
         set(hp1,'MarkerSize',6,'MarkerEdgeColor',[1 0 0], 'MarkerFacecolor', [1 0 0]);
         hp1 = plot(x2,0,'o');
         set(hp1,'MarkerSize',6,'MarkerEdgeColor',[1 0 0], 'MarkerFacecolor', [1 0 0]);
         xlabel('x','fontSize',12);
         ylabel('y','fontSize',12);
         title('Parabola','FontName','Arial','Fontsize',12);   
         hold off 
     end 
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  



