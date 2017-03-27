% cqify.m 
% Created: 03-23-2017 by JDR in Newark
% Last modified:
%
% Input: femStruct          - Strcut containing fem stuff
%        farFieldStruct     - Struct containing far field stuff
%        N                  - Number of centroids
%        MTime              - M+1 is the number of time steps
%        s                  - list of wavenumbers
%        t                  - time points  
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
% which is then transformed back. This code only works for time stepper
% methods such as BDF2 and Backward Euler.

function us = cqify(femStruct, farFieldStruct, N, MTime, s, t, lambda, flatP, P, ...
    iElements, jElements, nearFieldDistances, uiFun, qc, c0) 


%---- Generate scattered field data ----%
uiHat = sampleAndTransform(uiFun(femStruct.centroids(:,1),femStruct.centroids(:,2),t),lambda,MTime);

%-- Begin time-stepping routine. 
tic
uScatteredHat = timeStepper(N, MTime, s, farFieldStruct, femStruct,...
    P, flatP, iElements, jElements, uiHat, nearFieldDistances, qc, c0);
toc

% Calculate us in the time domain    
us = transformAndSample(uScatteredHat,lambda,MTime)-uiFun(femStruct.centroids(:,1),femStruct.centroids(:,2),t);



end