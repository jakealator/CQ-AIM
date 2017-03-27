% demoMoM.m 
% Created: 03-02-2017 by JDR in Newark
% Last Modified: 
%
% Runs a demo of the convolution quadrature-method of moments scheme
% to simulate acoustic wave propagation through an inhomogeneous medium.
% The geometry, numerical, and physical parameters are set in the files 
% called fowardParams.m and scatteringParams.m in ./demo. 

% To run this demo, simply click on Run or (if using Windows Shortcuts)
% type Shift+Enter

% For details about how this works and what it's doing, see README.md and 
% Chapter 5 of my Ph.D thesis.

%------ Begin Demo ------%


% Add required folders 
addpath(genpath('demo')) % Contains physical, geometric, and computational parameters for demo
addpath(genpath('modules')) % Contains the programs which actually compute

%---- Initialize parameters ----%
forwardParams
[femStruct, farFieldStruct, ~, ~, ~, ~, ~,flatP] = generateAuxillaryParams(meshStruct, N, M, d);
[iElements, jElements, multipoleMatrix] = generateNearFieldElements(N,0,0,0,0);


%---- Generate scattered field data ----%
X = femStruct.centroids;
nearFieldDistances = sqrt(bsxfun(@plus,full(dot(X',X',1)),full(dot(X',X',1))')-full(2*(X*X')));

uiHat = generateUIHat(uiHatFun, femStruct);


tic
qc = (1./(c(femStruct.centroids).^2)-1);
[KMat,MMat] = assembleNearFieldMatrices(femStruct.triAreas, nearFieldDistances, ...
        iElements, jElements, femStruct.centroids, zeros(N,N), qc,c0,waveNumber,N);

uHatMoM = (KMat+MMat)\(-(MMat)*(qc.*uiHat));
uScatteredHatMoM = uHatMoM - uiHat;

toc

%---- Plot results ----%
figure
subplot(2,1,1)
pltsln(meshStruct,femStruct.centroids,real(uScatteredHatMoM))
title('MoM - Real')

subplot(2,1,2)
pltsln(meshStruct,femStruct.centroids,imag(uScatteredHatMoM))
title('MoM - Imag')



%------ End Demo ------%