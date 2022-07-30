function DrawArrow(zT,magV,angleV,L,W,LW,col)

% Draw Arrow ==========================================================
% zT       position of vector tail (x,y): complex number zT = x + 1i*y
% magR     magnitude of vector
% angleV   angle of vector  [rad]
% L        length of arrow head
% W        width of arrow head
% LW       vector line width
% col       vector color

% (x,y) coordinates horizontal vector
  z1 = 0; z2 = magV;
  z3 = (magV-L) - 1i*W;
  z4 = conj(z3);
  Z0 = [z1 z2; z3 z4];
  
% Rotation matrix / rotation of vector
  rot = exp(1i*angleV);
  rotMatrix = [0 rot; rot rot];
  Z = Z0 .* rotMatrix;
  
% (x,y) position of vector tail
  xT = real(zT);
  yT = imag(zT);

% (x,y) poistion of vector head  
  x1 = real(Z(1,1)) + xT;
  y1 = imag(Z(1,1)) + yT;
  x2 = real(Z(1,2)) + xT;
  y2 = imag(Z(1,2)) + yT;
  x3 = real(Z(2,1)) + xT;
  y3 = imag(Z(2,1)) + yT;
  x4 = real(Z(2,2)) + xT;
  y4 = imag(Z(2,2)) + yT;

  X = [x1 x2 x2 x3 x2 x4];
  Y = [y1 y2 y2 y3 y2 y4];
  line(X,Y,'lineWidth',LW,'color',col)  
end

