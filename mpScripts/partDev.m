function pdf = partDev(a, x, w, eqType,n,m)

fPlus = zeros(n,m);
fMinus = zeros(n,m);
pdf = zeros(n,m);

for c = 1 : m
aP = a;
aM = a;
aP(c) = a(c)*1.001;
aM(c) = a(c)*0.999;

fPlus(:,c) = fitFunction(aP, x, eqType);
fMinus(:,c) = fitFunction(aM, x,eqType);
pdf(:,c) = (fPlus(:,c) - fMinus(:,c)) /( 0.002 * a(c));
end

for c = 1 : m                                   % weighted partial derivatives
pdf(:,c) = pdf(:,c) .* w';
end







