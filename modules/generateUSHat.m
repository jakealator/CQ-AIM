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
    jElements, extraFarFieldElements, qc, c0, Ng,N,waveNumber,P, flatGD)

% First compute near field components 
[KMat,MMat] = assembleNearFieldMatrices(triAreas, nearFieldDistances, ...
   iElements, jElements, extraFarFieldElements, qc, c0,waveNumber,N);
% KElements = -triAreas.^2.*(4*1i*c0^2/waveNumber^2+(2*pi*1i*c0/waveNumber).*...
%             besselh(1,1,1i*waveNumber/c0*sqrt(triAreas/pi)));
% KMat = sparse(iElements,jElements,KElements,N,N);
% KMat = (1i*waveNumber^2)/(4*c0^2)*KMat-waveNumber^2*extraFarFieldElements;
% MMat = sparse(1:N,1:N,triAreas./qc,N,N);

fftG = fft2(reshape(flatGD, 2*sqrt(Ng), 2*sqrt(Ng)));


% All we need to do now is compute the rhs and then use conjugate gradient to
% solve (I+V)x = rhs. 
%rhs = applyV((1./c(centroids).^2-1).*uiHat,KMat,fftG,P, waveNumber);
    rhs = MMat*uiHat;
    %usHat = gmres(@(x)applyIPlusV(x,MMat,KMat,fftG,P, waveNumber),rhs,30,1E-2);
    [usHat,~] = cgs(@(x)applyIPlusV(x,MMat,KMat,fftG,P, waveNumber),rhs,1E-3,20,@(x)(x./triAreas));
    


end

% aFarX applies operator Afar to a vector x. Output is approximation to Ax
% evaluated back on finite element mesh through the interpolation matrix. 
function Ax = applyIPlusV(x,M,K,fftG,P, waveNumber)

% The function apply V does the V application quickly, so we only need to
% add in the contribution from Mx
Ax = M*x+applyV(x,K,fftG, P, waveNumber);


end