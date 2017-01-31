% generateScatteredField.m
% Generates scattered data for time-dependent scattered from a penetrable
% inhomogeneous obstacle. Uses Galerkin-in-space Convolution
% Quadrature-in-time to generate data. 

% Created: 2015-04-21 by jdr at UD
% Last Modified: 

%function [uScatteredInt,uTotal,x,mesh] = generateScatteredField(transmitterLocations,... 
%    waveNumber, uIncident, c0, qc, mesh)

function [usScatteredExt, usScatteredInt] = generateScatteredFieldTest(s, c0, c, mesh, meshExt, d, uIncident)



% Initialize meshes (use p=1 even though that's not true)
mesh = initialize_mesh(mesh,1); % contains obstacles
meshExt = initialize_mesh(meshExt,1);
N=mesh.nt; % number of centroid points
NExt=meshExt.nt; % number of centroid points in exterior domain

qc=@(x,y)((c0^2./(c(x,y).^2)-1)+zeros(N,1)); % last part keeps size correct.

% Some values for computing on each mesh
triangleAreas = generateTriangleAreas(mesh,N);
[x,~,~,~] = generateCentroids(mesh,N);
[xExt,~,~,~] = generateCentroids(meshExt,NExt);

ui = zeros(N,1);
for k=1:N
    ui(k) = uIncident(d,x(k,:));
end
    
% Calculate matrices of differences
D = sqrt(bsxfun(@plus,full(dot(x',x',1)),full(dot(x',x',1))')-full(2*(x*x')));
DExt = sqrt(bsxfun(@plus,full(dot(x',x',1)),full(dot(xExt',xExt',1))')-full(2*(xExt*x')));

% Contrasts at centroid points
qj = qc(x(:,1),x(:,2));

H01 = besselh(0,1,1i*s/c0*D);
H11 = besselh(1,1,1i*s/c0*sqrt(triangleAreas/pi));
H01Ext = besselh(0,1,1i*s/c0*DExt);
    

[~,~,~, usScatteredInt] = generateUSHat(ui,triangleAreas,H01,H11,qj,c0,s,N,0);
%uTotal = uScatteredInt + ui;
usScatteredExt = generateUSHatExt(ui,usScatteredInt,triangleAreas,H01Ext,qj,s,c0,NExt);
end