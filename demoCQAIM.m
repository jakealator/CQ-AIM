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

%----Initialize parameters ----%
forwardParamsTime
[femStruct, farFieldStruct, iElements, jElements, multipoleMatrix, nearFieldDistances, P,flatP] = generateAuxillaryParams(meshStruct, N, M, d);
qc = (1./(c(femStruct.centroids).^2)-1);

%---- Generate scattered field data ----%

% Matrices containing data
uScatteredHat = zeros(N,M+1);
uiLambdaVec = zeros(N,M+1);

% Vectorize this!
for j=1:MTime+1
    uiLambdaVec(:,j) = lambda^(j-1)*uiFun(femStruct.centroids(:,1),femStruct.centroids(:,2),t(j));
end

uiHat = fft(uiLambdaVec,[],2);

%-- Begin time-stepping routine. Take advantage of symmetry of Fourier
%transform to half the number of needed solves
tic
for m=1:ceil((MTime-1)/2+1)
    
    m
    
    [flatGD, GD] = applyFundamentalSolution(s(m), farFieldStruct.farFieldGrid);

    extraFarFieldElements = assembleFarFieldMatrix(s(m),  flatP, N, ...
    farFieldStruct.rectangularElementsX, farFieldStruct.rectangularElementsY, ...
    iElements, jElements, farFieldStruct.rectangularLocations, GD);

    uScatteredHat(:,m) = generateUSHat(uiHat(:,m), femStruct.triAreas, nearFieldDistances, iElements, ...
    jElements, femStruct.centroids, extraFarFieldElements, c, c0, farFieldStruct.nG,N,s(m), P, flatGD, farFieldStruct);
end

uScatteredHat = 1./(2*qc).*uScatteredHat;
uScatteredHat(:,(ceil((MTime-1)/2+1)+1):(MTime+1)) = conj(uScatteredHat(:,ceil((MTime-1)/2+1):-1:2));

% Calculate us in the time domain    
us = ifft(uScatteredHat,[],2);
for m=1:MTime+1
    us(:,m) = (lambda^(-m+1))*us(:,m); 
end
us = real(us); 

toc

%---- Plot results ----%
figure
for j=1:201
pltsln(meshStruct,femStruct.centroids,us(:,j))
%axis([-0.5,0.5,-0.5,0.5,min(min(us)),max(max(us))])
caxis([min(min(us(:,1:100))),max(max(us(:,1:100)))])
title(t(j))
view(2)
pause(0.1)

end


%------ End Demo ------%




