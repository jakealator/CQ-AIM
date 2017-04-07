

function uScatteredHat = timeStepper(N, MTime, s, farFieldStruct, femStruct,...
    P, flatP, iElements, jElements, uiHat, nearFieldDistances, qc, c0) 

uScatteredHat = zeros(N,MTime+1);

% inds = 1:ceil((MTime-1)/2+1);
% solveInd = inds(~(max(abs(uiHat(:,inds)))<1E-2));

numberSolves = 0;
parfor m=1:ceil((MTime-1)/2+1)
% parfor m = solveInd


            if max(abs(uiHat(:,m)))<1E-2
                uScatteredHat(:,m) = zeros(N,1);
            else

                numberSolves = numberSolves+1;
                [flatGD, GD] = applyFundamentalSolution(s(m), farFieldStruct.farFieldGrid);

                extraFarFieldElements = assembleFarFieldMatrix(s(m),  flatP, N, ...
                farFieldStruct.rectangularElementsX, farFieldStruct.rectangularElementsY, ...
                iElements, jElements, farFieldStruct.rectangularLocations, GD);

                uScatteredHat(:,m) = generateUSHat(uiHat(:,m), femStruct.triAreas, nearFieldDistances, iElements, ...
                jElements, extraFarFieldElements, qc, c0, farFieldStruct.nG,N,s(m), P, flatGD);
            end


end

% uScatteredHat = uScatteredHat;
uScatteredHat(:,(ceil((MTime-1)/2+1)+1):(MTime+1)) = conj(uScatteredHat(:,ceil((MTime-1)/2+1):-1:2));

sprintf(strjoin({'Time stepping complete after ', int2str(numberSolves), ' solves.'}))


end

% For full mom

%         [iElements, jElements, ~] = generateNearFieldElements(N,0,0,0,0);
%         X = femStruct.centroids;
%         nearFieldDistances = sqrt(bsxfun(@plus,full(dot(X',X',1)),full(dot(X',X',1))')-full(2*(X*X')));
%         [KMat,MMat] = assembleNearFieldMatrices(femStruct.triAreas, nearFieldDistances, ...
%             iElements, jElements, zeros(N,N), qc,c0,s(m),N);
%         uScatteredHat(:,m) = (KMat+MMat)\((MMat)*(uiHat(:,m)));
