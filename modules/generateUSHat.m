% generateUSHat.m
% Created on: 02-17-2017 by JDR in Newark
% Last Modified: 02-20-2017
%
% Inputs: uiHat                 - the Nx1 incident field corresponding to current
%                                 wavenumber 
%         triAreas              - Nx1 vector of triangle areas
%         nearFieldDistances    - Distance matrix for near field distances
%         i/jElements           - Vectors indicating which elements of FEM
%                                 mesh are near field 
%         centroids             - Nx2 vector of centroids
%         extraFarFieldElements - Far field components which need to be
%                                 subtracted away (same size as i/jElements)
%         c/c0                  - speed of sound inside/ouside medium
%         Ng                    - Number of Cartesian grid points
%         N                     - Number of FEM grid points
%         waveNumber            - current wavenumber
%         P                     - NxNg interpolation matrix
%
% Outputs: usHat - Nx1 vector containing the uSHat data for a fixed
%                  wavenumber. 
%
% Generates the scattered field corresponding to a wavenumber k, an index
% of refraction defined by c and the shape D, and an incident field ui.
% This is just a time harmonic scattered field. However, it is computed in
% the near field with a P0-Galerkin approximation to the Lippmann-Schwinger
% equation and in the far field with a FFT-based techinque, using
% delta-function 'basis functions' defined on a rectangular grid. 
%
% Algorithm: 1) Calculate near field based on previously-computed near
%               field grid. 
%            2) Compute RHS
%            3) Use CG to solve for US in D, using 1) and far field
%               computed during iterations
%            4) Find solution in exterior. 

function usHat = generateUSHat(uiHat,triAreas, nearFieldDistances, iElements, ...
    jElements, centroids, extraFarFieldElements, c, c0, Ng,N,waveNumber,flatP,P, GD, farFieldStruct)

% First compute near field components 
[KMat,MMat] = assembleNearFieldMatrices(triAreas, nearFieldDistances, ...
    iElements, jElements, centroids, extraFarFieldElements, c,c0,waveNumber,N);

% Add in fft of Green's function at rectangular points
% GdLong = [Gd; zeros(length(Gd),1)];
% fftG = fft(GdLong);
% % zero pad
GdLong = [zeros(length(GD), 2*length(GD)); GD, zeros(length(GD))];
fftG = fft2(GdLong);
%  fftG = generateFullFarFieldMatrix(zeros(N,1),farFieldStruct,flatP, waveNumber);



% All we need to do now is compute the rhs and then use conjugate gradient to
% solve (I+V)x = rhs. 
rhs = -waveNumber*applyV((1./c(centroids).^2-1).*uiHat,KMat,N,fftG,P, waveNumber, c0, farFieldStruct);
usHat = cgs(@(x)applyIPlusV(x,MMat,KMat,N,fftG,P, waveNumber, c0, farFieldStruct),rhs,1E-4);


end

% aFarX applies operator Afar to a vector x. Output is approximation to Ax
% evaluated back on finite element mesh through the interpolation matrix. 
function Ax = applyIPlusV(x,M,K,N,fftG,P, waveNumber, c0, farFieldStruct)

% The function apply V does the V application quickly, so we only need to
% add in the contribution from Mx
Ax = M*x+applyV(x,K,N,fftG, P, waveNumber, c0, farFieldStruct);



end