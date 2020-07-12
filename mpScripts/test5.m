clf
N=50;
[X,Y,Z]=sphere(N);
C=ones(size(X));
az=160;el=20;
%%%% Plot the sphere with coordinate lines %%%%%%%%%%
rjxsurf(X,Y,Z,C,'setlines','on','az',az,'el',el,'linesep1',2,...
'linesep2',2,'col',[0 0 1],'linecol',[1 1 0]);
%%%% Set out a ’surface line’ on the sphere %%%%%%%%%%
phibase=linspace(0,2*pi,N+1); % A coordinate basevector
thetabase=linspace(0,pi,N+1); % A coordinate basevector
phi=linspace(2.8,5.2,200); % Trajectory phi
theta=pi/2+sin(phi*10).*(phi-min(phi))*0.2; % Trajectory theta
rjcell{1}={'surfline',theta,phi,'col',[1 0 1],'epsilon',0.1,...
'arrow','curved'};
handles = putonsurf(X,Y,Z,thetabase,phibase,{az,el,[]},rjcell);
rjsurfset(handles,{0.7,0.3,0.4})
%%%% Fixing some lighting etc %%%%%%%%%%
view(az,el);
axis off;axis image
shading flat
light;lighting phong