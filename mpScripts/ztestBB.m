
S = spinop2('gl')

dom = 50*[-1 1 -1 1];
tspan = [0 16];
S = spinop2(dom,tspan);
S.lin = @(u) lap(u);
S.nonlin = @(u) u - (1+1.5i)*u.*(abs(u).^2);

x = chebfun2(@(x,y) x,dom); y = chebfun2(@(x,y) y,dom);
u1 = (1i*x+y).*exp(-.03*(x.^2+y.^2)); S.init = u1;
npts = 80; dt = 4/npts; tic
u = spin2(S,npts,dt,'plot','off');
plot(real(u)), view(0,90), axis equal, axis off