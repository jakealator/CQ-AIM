% cqify.m 
% Created: 03-23-2017 by JDR in Newark
% Last modified:
%
% Input: femStruct          - Strcut containing fem stuff
%        farFieldStruct     - Struct containing far field stuff
%        N                  - Number of centroids
%        MTime              - M+1 is the number of time steps
%        t                  - time points 
%        dt                 - time stpe
%        lambda             - parameter related to accuracy of approximate
%                             CQ
%        flatP              - nGx9 matrix containing nonzero elements
%                             mapping from fem grid to far field grid 
%        P                  - NxnG matrix of the expanded version of above
%        i/jElements        - Contain locations of near field elements   
%        nearFieldDistances - Distances between near field elements
%        uiFun              - Function to compute incident field
%        c                  - Function to commpute speed of sound in medium
%        c0                 - Speed of sound in free space (needs to be 1)
%        qc                 - Contrast, vector of size Nx1
% Output: us - Nx(M+1) vector of scattered field
%
% Uses convolution quadrature to compute scattered field using
% timeStepper.m. In particular, the timeStepper.m code uses the algorithms
% from chapter 4 of Hassel-Sayas along with a fast AIM Helmholtz Equation
% solver to compute the Fourier domain version of the scattered field,
% which is then transformed back. This code only works for RK time step
% methods. 

function us = RKcqify(femStruct, farFieldStruct, N, MTime, t, dt, lambda, flatP, P, ...
    iElements, jElements, nearFieldDistances, uiFun, qc, c0) 

A=[5/12 -1/12; 3/4 1/4]; cRK = [1/3,1]; % Radau-II
%A = [1/6 -1/3 1/6; 1/6 5/12 -1/12; 1/6 2/3 1/6]; cRK = [0 1/2 1]; % Radau-III
%---- Generate scattered field data ----%
    numberStages = length(A);
    uiSample = zeros(N*numberStages,MTime+1);
    for p=1:numberStages
        currentElements = (p-1)*N+1:p*N;
        uiSample(currentElements,:) = uiFun(femStruct.centroids(:,1),femStruct.centroids(:,2),t+dt*cRK(p));
    end
 uiHat = sampleAndTransform(uiSample,lambda,MTime);

%-- Begin time-stepping routine. 
tic
uScatteredHat = RKtimeStepper(N, MTime, A, farFieldStruct, femStruct,...
    P, flatP, iElements, jElements, uiHat, nearFieldDistances, qc, c0, lambda,dt);
toc

% Calculate us in the time domain  
usFull = transformAndSample(uScatteredHat,lambda,MTime);
us = usFull((length(A)-1)*N+1:end,:);


end