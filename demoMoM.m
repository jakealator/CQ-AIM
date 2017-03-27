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
forwardParamsTime
[femStruct, farFieldStruct, ~, ~, ~, ~, ~,flatP] = generateAuxillaryParams(meshStruct, N, M, d);
[iElements, jElements, multipoleMatrix] = generateNearFieldElements(N,0,0,0,0);


%---- Generate scattered field data ----%
X = femStruct.centroids;
nearFieldDistances = sqrt(bsxfun(@plus,full(dot(X',X',1)),full(dot(X',X',1))')-full(2*(X*X')));

% extraFarFieldElements = waveNumber^2*assembleFarFieldMatrix(waveNumber,  flatP, N, ...
%     farFieldStruct.rectangularElementsX, farFieldStruct.rectangularElementsY, ...
%     iElements, jElements);

tic

% Matrices containing data
uScatteredHatMoM = zeros(N,M+1);
uiLambdaVec = zeros(N,M+1);

% Vectorize this!
for j=1:MTime+1
    uiLambdaVec(:,j) = lambda^(j-1)*uiFun(femStruct.centroids(:,1),femStruct.centroids(:,2),t(j));
end

uiHat = fft(uiLambdaVec,[],2);

%-- Begin time-stepping routine. Take advantage of symmetry of Fourier
%transform to half the number of needed solves
for m=1:ceil((MTime-1)/2+1)
    if ~mod(m,10)
        m
    end
    
    [KMat,MMat] = assembleNearFieldMatrices(femStruct.triAreas, nearFieldDistances, ...
        iElements, jElements, femStruct.centroids, zeros(N,N), c,c0,s(m),N);


    %uScatteredHatMoM(:,m) = (KMat+MMat)\(-(KMat)*((1./c(femStruct.centroids).^2-1).*uiHat(:,m)));
    uScatteredHatMoM(:,m) = gmres((KMat+MMat),-(KMat)*((1./c(femStruct.centroids).^2-1).*uiHat(:,m)),10,1E-3,50);
end
uScatteredHatMoM = 1./((1./c(femStruct.centroids).^2-1)).*uScatteredHatMoM;
uScatteredHatMoM(:,(ceil((MTime-1)/2+1)+1):(MTime+1)) = conj(uScatteredHatMoM(:,ceil((MTime-1)/2+1):-1:2));

% Calculate us in the time domain    
usMoM = ifft(uScatteredHatMoM,[],2);
for m=1:MTime+1
    usMoM(:,m) = (lambda^(-m+1))*usMoM(:,m); 
end
usMoM = real(usMoM); % NOTE THIS!

toc

%---- Plot results ----%
figure
for j=1:MTime+1
pltsln(meshStruct,femStruct.centroids,usMoM(:,j))
%axis([-0.5,0.5,-0.5,0.5,min(min(us)),max(max(us))])
caxis([min(min(usMoM)),max(max(usMoM))])
title(t(j))
view(2)
pause(0.1)

end






%------ End Demo ------%