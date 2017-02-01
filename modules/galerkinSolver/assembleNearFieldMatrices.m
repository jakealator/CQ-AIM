% assembleNearFieldMatrices.m
% Created: 02-01-2017 by JDR at UD
% Last Modified: 
%
% Input:
%
% Output:
%
% Outputs the near field components of the AIM algorithm for computing the
% action of the time-harmonic Lippmann-Schwinger equation on a vector. Uses
% p=0 piecewise constant elements on each element and approximates
% integrals on triangles with an exact formula on an equal area circle. 

function [K,M,nearFieldMatrix] = assembleNearFieldMatrices(triAreas,nearFieldDistances, iElements, jElements, centroids, c,c0,s,N)

% Contrasts at centroid points
qj =((c0./c(centroids(:,1),centroids(:,2))).^2-ones(N,1));


KElements = zeros(length(iElements),1); % Will hold values of K ('stiffness') matrix
MElements = zeros(length(iElements),1); % Will hold values of M (mass) matrix
%K(diagonal) = -triangleAreas.*(4*1i*c0^2/s^2+(2*pi*1i*c0/s).*H11);
%K(off-diagonal) = triangleAreas(i)*triangleAreas(j)*H01(i,j)
% Only calculate values where we have specified elements to be 'near field'
for k = 1:length(iElements)
    if iElements(k)==jElements(k) % On the diagonal
        KElements(k) = -triAreas(iElements(k))^2.*(4*1i*c0^2/s^2+(2*pi*1i*c0/s).*...
            besselh(1,1,1i*s/c0*sqrt(triAreas(iElements(k))/pi)));
        MElements(k) = triAreas(iElements(k))/qj(iElements(k));
    else
        KElements(k) = triAreas(iElements(k))*triAreas(jElements(k))*...
            besselh(0,1,1i*s/c0*nearFieldDistances(iElements(k),jElements(k)));
    end
end
    
% Create a sparse matrix whose diagonal elements correspond to qj at
% required are exactly the ones
% created above 
M = sparse(iElements,jElements,MElements,N,N);
K = sparse(iElements,jElements,KElements,N,N);

K = (1i*s^2)/(4*c0^2)*K;

nearFieldMatrix = M+K;

end