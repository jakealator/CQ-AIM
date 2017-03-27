

function uScatteredHat = timeStepper(N, MTime, s, farFieldStruct, femStruct,...
    P, flatP, iElements, jElements, uiHat, nearFieldDistances, qc, c0) 

uScatteredHat = zeros(N,MTime+1);

parfor m=1:ceil((MTime-1)/2+1)
    
    m
    abs(s(m))
%     if abs(s(m))<0
%         0
%         flatP = generateInterpolationMatrix(farFieldStruct.centers,...
%             farFieldStruct.rectangularElementsX, farFieldStruct.rectangularElementsY,...
%             0,0, 0, N,s(m),femStruct.centroids);
% 
%         femIndex = repmat(1:N,9,1);
%         femIndex = femIndex(:);
%         rectangularLocationsSideways=farFieldStruct.rectangularLocations.';
%         rectangularLocationsSideways = rectangularLocationsSideways(:);
%         flatPSideways = flatP.';
%         P = sparse(rectangularLocationsSideways,femIndex,flatPSideways(:),farFieldStruct.nG,N); 
%     end
%         
        [flatGD, GD] = applyFundamentalSolution(s(m), farFieldStruct.farFieldGrid);

        extraFarFieldElements = assembleFarFieldMatrix(s(m),  flatP, N, ...
        farFieldStruct.rectangularElementsX, farFieldStruct.rectangularElementsY, ...
        iElements, jElements, farFieldStruct.rectangularLocations, GD);

        uScatteredHat(:,m) = generateUSHat(uiHat(:,m), femStruct.triAreas, nearFieldDistances, iElements, ...
        jElements, extraFarFieldElements, qc, c0, farFieldStruct.nG,N,s(m), P, flatGD);
    
%         1
%         [iElements, jElements, ~] = generateNearFieldElements(N,0,0,0,0);
%         X = femStruct.centroids;
%         nearFieldDistances = sqrt(bsxfun(@plus,full(dot(X',X',1)),full(dot(X',X',1))')-full(2*(X*X')));
%         [KMat,MMat] = assembleNearFieldMatrices(femStruct.triAreas, nearFieldDistances, ...
%             iElements, jElements, zeros(N,N), qc,c0,s(m),N);
%         uScatteredHat(:,m) = (KMat+MMat)\((MMat)*(uiHat(:,m)));

end

% uScatteredHat = uScatteredHat;
uScatteredHat(:,(ceil((MTime-1)/2+1)+1):(MTime+1)) = conj(uScatteredHat(:,ceil((MTime-1)/2+1):-1:2));




end