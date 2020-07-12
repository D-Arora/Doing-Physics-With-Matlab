clc; close all; clear; a=input('Enter the sphere radius a <= 2 [m] = ');
E0=input('Enter the intensity of external uniform E-field [V/m]: ');
Eps0=8.854187817e-12; Phi=linspace(0,pi,18); N=2e2; D=6; dd=linspace(-D/2,D/2,N);
[X,Y,Z]=meshgrid(dd); [az,el,r] = cart2sph(X,Y,Z); U=-E0*(1-a^3./r.^3).*r.*sin(el); Ud=E0*a^3./r.^3.*r.*sin(el);
[Ex,Ey,Ez]=gradient(-U,D/N); E=sqrt(Ex.^2+Ey.^2+Ez.^2); FF=20*log10(abs(U));
[Exd,Eyd,Ezd]=gradient(-Ud,D/N); Ed=sqrt(Exd.^2+Ey.^2+Ez.^2); FFd=Ud;
f1=figure('units','normalized','outerposition',[0 0 1 1]);
hold on; grid on; view(-138,9); axis square; A=[-1 1 -1 1 -1 1]*D/2; axis(A);
[xs,ys,zs] = sphere(50); xs=xs*a; ys=ys*a; zs=zs*a; hC=surf(xs,ys,zs,'FaceColor','w');
r0 = a; b = 1 ; h = a/10; X1 = xs ; Y1 = ys ; Z1 = zs ; X1(zs<r0/2) = NaN ; Y1(zs<r0/2) = NaN ; Z1(zs<r0/2) = NaN ;
surf(X1,Y1,Z1,'FaceColor','r'); X2 = xs ; Y2 = ys ; Z2 = zs ; X2(zs>-r0/2) = NaN ; Y2(zs>-r0/2) = NaN ; Z2(zs>-r0/2) = NaN ;
surf(X2,Y2,Z2,'FaceColor','b'); h1=slice(X,Y,Z,FF,[],0,[]); set(h1,'FaceColor','interp','EdgeColor','none');
hs1=streamslice(X,Y,Z,Ex,Ey,Ez,[],0,[],2); set(hs1,'Color','r');
c=hsv; c=c(1:2:end,:); colormap(c); hc=colorbar; ylabel(hc,'\bf 20*log10(Electric Potential Intensity) [dB]','FontSize',20)
zlabel('\bfZ-axis'); ylabel('\bfY-axis'); xlabel('\bfX-axis'); zArrow = [0-D/2 1 1 1.2 1 1 0-D/2]*(D/2)*0.8;
yArrow = [-0.1 -0.1 -0.2 0 0.2 0.1 0.1]*D/2; xArrow = 0*yArrow; hArrow = fill3(-D/1e3+xArrow,yArrow+3*D/8,zArrow,[1 0 0]);
hArrow = fill3(-D/1e3+xArrow,yArrow-3*D/8,zArrow,[1 0 0]); hArrow = fill3(xArrow+3*D/8,yArrow+3*D/8,zArrow,[1 0 0]);
hArrow = fill3(xArrow+3*D/8,yArrow-3*D/8,zArrow,[1 0 0]); for jj=2:length(Phi); dT=(jj-1)*(Phi(jj)-Phi(jj-1))*180/pi;
h1=slice(X,Y,Z,FF,0,[],[]); set(h1,'FaceColor','interp','EdgeColor','none');
hs1=streamslice(X,Y,Z,Ex,Ey,Ez,0,[],[],2); set(hs1,'Color','r'); rotate(h1,[0,0,1],dT);
rotate(hs1,[0,0,1],dT); drawnow; set(h1,'visible','off'); set(hs1,'visible','off'); end
phi=0:pi/100:pi; phi1=phi/4; [U1,V1]=pol2cart(phi,1+sin(phi)); [U2,V2]=pol2cart(-phi,1+sin(phi));
[U3,V3]=pol2cart(phi,1+0*phi); [U4,V4]=pol2cart(-phi,1+0*phi); [U5,V5]=pol2cart(phi1,0.5+0*phi1);
f2=figure; grid minor; axis equal;movegui(f2,'west'); hold on; grid on; view(-138,9); axis square; A=[-1 1 -1 1 -1 1]*D/2; axis(A);
[xs,ys,zs] = sphere(50); xs=xs*a; ys=ys*a; zs=zs*a; hC=surf(xs,ys,zs,'FaceColor','w');
r0 = a; b = 1 ; h = a/10; X1 = xs ; Y1 = ys ; Z1 = zs ; X1(zs<r0/2) = NaN ; Y1(zs<r0/2) = NaN ; Z1(zs<r0/2) = NaN ;
surf(X1,Y1,Z1,'FaceColor','r'); X2 = xs ; Y2 = ys ; Z2 = zs ; X2(zs>-r0/2) = NaN ; Y2(zs>-r0/2) = NaN ; Z2(zs>-r0/2) = NaN ;
surf(X2,Y2,Z2,'FaceColor','b'); h1=slice(X,Y,Z,FFd,[],0,[]); set(h1,'FaceColor','interp','EdgeColor','none');
hs1=streamslice(X,Y,Z,Exd,Eyd,Ezd,[],0,[],2); set(hs1,'Color','r');
c=hsv; c=c(1:2:end,:); colormap(c); hc=colorbar; ylabel(hc,'\bf Electric Potential Intensity [V]','FontSize',20)
zlabel('\bfZ-axis'); ylabel('\bfY-axis'); xlabel('\bfX-axis'); title('\bfSelf Potential and E-force Lines of Sphere'); movegui(f2,'west');
f3=figure; grid minor; axis equal; movegui(f3,'east'); [xs,ys,zs] = sphere(50); xs=xs*a; ys=ys*a; zs=zs*a;
hC=surf(xs,ys,zs);[azimuth,elevation,r1] = cart2sph(xs,ys,zs);
Sigma=sin(elevation); set(hC,'CData',Sigma); colormap('jet'); shading interp;
title('\bfNormilized Surface Charge on Sphere'); colorbar