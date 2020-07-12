% pc_001.m

   clear all
   close all
   clc
   
% Inputs
   L = 32;
   p = 0.50;
% Setup
   RN = rand(L,L);
   u = zeros(L,L);
   u(RN < p) =  1;
 
% clusters:  cls = cluster labels 1, 2, ...   numC = number of clusters 
   [cls, numC] = bwlabel(u,4);
% number of sites
   nSites  = L*L; 
% number of occupied sites
   nOSites = sum(sum(u == 1));
% probability of a site being occupied
   pOSites  = nOSites/nSites;     

% number of occupied sites in each cluster   
   nC = zeros(numC,1);  
   for c = 1 : numC
     nC(c) = length(cls(cls == c)); 
   end
% s   counter from 1 to max number of occupied sites in a cluster  
    s = 1:max(nC);
% nS  cluster size distribution 
    nS = zeros(max(nC),1);
  for c = s
   nS(c) = length(nC(nC == c));
  end
% Probability that an occupied site chosen at random is part of an s-site cluster
   probROS = nS ./ sum(nS);
% mean clustert size = mean number of occupied sites in cluster
   meanClusterSize = mean(nC);    
   S = sum(s.*probROS');
   
   
% CHECKING FOR SPANNING CLUSTERS
%    X spanning
   isSpanX = zeros(numC,1); isSpanY = isSpanX;
   %isSpanX = zeros(numC,1); isSpanY = isSpanX;
for c = 1 : numC
    col1 = isempty(find(cls(cls(:,1)==c)));
    colN = isempty(find(cls(cls(:,end)==c)));
    col1N = col1 + colN;
    if col1N == 0; isSpanX(c) = 1;end
    col1 = isempty(find(cls(cls(1,:)==c)));
    colN = isempty(find(cls(cls(end,:)==c)));
    col1N = col1 + colN;
    if col1N == 0; isSpanY(c) = 1;end
    
end

   
% OUTPUT   ==============================================================
disp('  ');
fprintf('Linear dimension  L = %2.0f \n',L);
fprintf('Total number of sites  nSties = %2.0f \n',nSites);
fprintf('Probability of a site being occupied   p = %2.2f \n',p);
disp('  ');
fprintf('Number of occupied sites nOSites = %2.0f \n',nOSites);
fprintf('Simulation: probability of a site being occupied pOSites = %2.2f \n',pOSites);
fprintf('Number of clusters  numC = %2.0f \n',numC);
disp('Cluster Size: number of occupied sites in each cluster  nC ');
fprintf('   %2.0f ',nC');
disp('  ')
disp('Cluster Size Distribution  ');
disp('    s      nS    probROS');
for c = s
 fprintf('   %2.0f     %2.0f     %2.3f   \n',s(c), nS(c), probROS(c));  
end
disp('  ');
disp('Mean Cluster Size  S:  sites / clusters / prob   ')
fprintf('S = mean(nC) = %2.3f \n',meanClusterSize);
fprintf('S = SUM(s .* probROS) = %2.3f \n',S);
disp('  ');
disp('Spanning clusters = 1')
disp('cluster #   rows   cols'); 

for c = 1: numC
  fprintf('     %2.0f      %2.0f       %2.0f     \n',c, isSpanY(c), isSpanX(c));  
end
   
 
%   GRAPHICS =============================================================
   figure(1)
   set(gcf,'units','normalized','position',[0.1 0.52 0.23 0.32]);
   hold on
   for cx = 1 : L
       for cy = 1 : L
       if u(cy,cx) == 0
       h_rectangle = rectangle('Position',[cx-0.5,cy-0.5,1,1]);    
       set(h_rectangle,'faceColor',[1 1 1]);
       end
              
       if u(cy,cx) == 1
       h_rectangle = rectangle('Position',[cx-0.5,cy-0.5,1,1]);    
       set(h_rectangle,'faceColor',[0 0 1],'edgeColor',[0 0 1]);
       end
       
       
       end
   end
   
   for c = 1 : numC;
      [rN, cN] = find(cls == c);
      text(cN(1),rN(1),num2str(c,'%2.0f'),'Color','w');
   end
   
   axis equal
   axis tight
   
% Color coding clusters ==================================================
   clsT = zeros(L,L);
   for c = 1 : L
      clsT(c,:) = cls(L+1-c,:);
   end

hf2 = figure(2);
   set(gcf,'units','normalized','position',[0.35 0.52 0.23 0.32]);
    
   img = label2rgb(clsT);
   image(img);
    hold on
    for cx = 0 : L
     for cy = 0 : L 
         yP = 0.5:L+0.5; xP = 0.5 + cx .* ones(length(yP),1);
         plot(xP,yP,'k','linewidth',0.5);
         xP = 0.5:L+0.5; yP = 0.5 + cx .* ones(length(yP),1);
         plot(xP,yP,'k','linewidth',0.5);
     end
    end
    axis off
    box on
%    image(img)
  