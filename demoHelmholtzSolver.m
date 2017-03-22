% demoHelmholtzSolver.m 
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

%----Initialize parameters ----%
forwardParams
[femStruct, farFieldStruct, iElements, jElements, multipoleMatrix, nearFieldDistances, P,flatP] = generateAuxillaryParams(meshStruct, N, M, d);
uiHat = generateUIHat(uiHatFun, femStruct);

%---- Generate scattered field data ----%
tic

[flatGD, GD] = applyFundamentalSolution(waveNumber, farFieldStruct.farFieldGrid);

extraFarFieldElements = ...
 assembleFarFieldMatrix(waveNumber,  flatP, N, ...
     farFieldStruct.rectangularElementsX, farFieldStruct.rectangularElementsY, ...
     iElements, jElements, farFieldStruct.rectangularLocations, GD);

uScatteredHat = generateUSHat(uiHat, femStruct.triAreas, nearFieldDistances, iElements, ...
jElements, femStruct.centroids, extraFarFieldElements, c, c0, farFieldStruct.nG,N,waveNumber, P, flatGD, farFieldStruct,zeros(N,1));
qc = (1./(c(femStruct.centroids).^2)-1); uScatteredHat = 1./qc*uScatteredHat;


toc


%---- Plot results ----%
uScatteredHat = 1/2*uScatteredHat;
figure
subplot(2,1,1)
pltsln(meshStruct,femStruct.centroids,real(uScatteredHat))
title('AIM - Real')

subplot(2,1,2)
pltsln(meshStruct,femStruct.centroids,imag(uScatteredHat))
title('AIM - Imag')

%------ End Demo ------%