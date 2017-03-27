



function uScatteredHat = RKtimeStepper(N, MTime, A, farFieldStruct, femStruct,...
    P, flatP, iElements, jElements, uiHat, nearFieldDistances, qc, c0, lambda,dt) 


numberStages = length(A);
% For stiff RK methods, delta = A^{-1}-zeta*A^{-1}1.e_{numberStages}
Ainv = inv(A);
C = Ainv*ones(length(A),1)*[zeros(numberStages-1,1)',1];

xi = exp((-2*pi*1i*(0:MTime))/(MTime+1));

% uScatteredHat contains all of the stages of the RK process - for the mth
% time step, uScatteredHat((numberSteps-1)*N+1:numberSteps*N,m) are the
% stages, with the last numberSteps containing the actual time step. 
uScatteredHat = zeros(N*numberStages,MTime+1);
for m=1:ceil((MTime-1)/2+1)
    m
    [PEig,Lambda]=eig(1/dt*(Ainv-lambda*xi(m)*C));
    Lambda=diag(Lambda);
    rhs = kron(inv(PEig),speye(N))*uiHat(:,m);
    for p=1:numberStages
        currentElements = (p-1)*N+1:p*N;
        [flatGD, GD] = applyFundamentalSolution(Lambda(p), farFieldStruct.farFieldGrid);

        extraFarFieldElements = assembleFarFieldMatrix(Lambda(p),  flatP, N, ...
            farFieldStruct.rectangularElementsX, farFieldStruct.rectangularElementsY, ...
            iElements, jElements, farFieldStruct.rectangularLocations, GD);

        uScatteredHat(currentElements, m) = generateUSHat(rhs(currentElements), femStruct.triAreas, nearFieldDistances, iElements, ...
            jElements, extraFarFieldElements, qc, c0, farFieldStruct.nG,N,Lambda(p), P, flatGD);
        
    end
    uScatteredHat(:,m)=kron(PEig,speye(N))*uScatteredHat(:,m);
    
end

% uScatteredHat = 1./(2*qc.^2).*uScatteredHat;
uScatteredHat(:,(ceil((MTime-1)/2+1)+1):(MTime+1)) = conj(uScatteredHat(:,ceil((MTime-1)/2+1):-1:2));




end