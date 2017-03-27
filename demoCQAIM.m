% demoCQAIM.m 
% Created: 01-26-2016 by JDR in Newark
%
% Runs a demo of the convolution quadrature-adaptive integral method scheme
% to simulate acoustic wave propagation through an inhomogeneous medium.
% The geometry, numerical, and physical parameters are set in the files 
% called fowardParams.m and scatteringParams.m in ./demo. 

% To run this demo, simply click on Run or (if using Windows Shortcuts)
% type Shift+Enter

% For details about how this works and what it's doing, see README.md and 
% Chapter 5 of my Ph.D thesis.

%------ Begin Demo ------%
clear

% Add required folders 
addpath(genpath('demo')) % Contains physical, geometric, and computational parameters for demo
addpath(genpath('modules')) % Contains the programs which actually compute

%---- Initialize parameters ----%
forwardParamsTime
[femStruct, farFieldStruct, iElements, jElements, multipoleMatrix, nearFieldDistances, P,flatP] = generateAuxillaryParams(meshStruct, N, M, d);
qc = (1./(c(femStruct.centroids).^2)-1);

tic
%---- Find scattered Field ---%
%us = RKcqify(femStruct, farFieldStruct, N, MTime, t, dt, lambda, flatP, P, ...
%    iElements, jElements, nearFieldDistances, uiFun, c, c0, qc);
 us = cqify(femStruct, farFieldStruct, N, MTime, s, t, lambda, flatP, P, ...
     iElements, jElements, nearFieldDistances, uiFun, qc, c0);

toc

%---- Plot results ----%
figure
for j=1:MTime+1
pltsln(meshStruct,femStruct.centroids,us(:,j))
axis([-0.5,0.5,-0.5,0.5,min(min(us)),max(max(us))])
 caxis([min(min(us)),max(max(us))])
title(t(j))
% view(2)
pause(0.1)

end


%------ End Demo ------%




