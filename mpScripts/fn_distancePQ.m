function d = distancePQ(xP,yP,zP,xQ,yQ,zQ)

% Function to calculate the distance between two points P and Q

 d = sqrt((xP-xQ).^2 + (yP-yQ).^2 + (zP-zQ)^2);
end

% function d = distance(xP,yP,zP,xQ,yQ,zQ)
% 

% 
% d = sqrt((xP-xQ)^2 + (yP-yQ)^2 + (zP-zQ)^2);